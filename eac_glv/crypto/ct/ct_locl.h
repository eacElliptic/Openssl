/*
 * Written by Rob Percival (robpercival@google.com) for the OpenSSL project.
 */
/* ====================================================================
 * Copyright (c) 2016 The OpenSSL Project.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. All advertising materials mentioning features or use of this
 *    software must display the following acknowledgment:
 *    "This product includes software developed by the OpenSSL Project
 *    for use in the OpenSSL Toolkit. (http://www.OpenSSL.org/)"
 *
 * 4. The names "OpenSSL Toolkit" and "OpenSSL Project" must not be used to
 *    endorse or promote products derived from this software without
 *    prior written permission. For written permission, please contact
 *    licensing@OpenSSL.org.
 *
 * 5. Products derived from this software may not be called "OpenSSL"
 *    nor may "OpenSSL" appear in their names without prior written
 *    permission of the OpenSSL Project.
 *
 * 6. Redistributions of any form whatsoever must retain the following
 *    acknowledgment:
 *    "This product includes software developed by the OpenSSL Project
 *    for use in the OpenSSL Toolkit (http://www.OpenSSL.org/)"
 *
 * THIS SOFTWARE IS PROVIDED BY THE OpenSSL PROJECT ``AS IS'' AND ANY
 * EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE OpenSSL PROJECT OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 * ====================================================================
 */

#ifdef OPENSSL_NO_CT
# error CT is disabled.
#endif

#include <stddef.h>

#include <openssl/ct.h>
#include <openssl/evp.h>
#include <openssl/x509.h>
#include <openssl/safestack.h>

/*
 * From RFC6962: opaque SerializedSCT<1..2^16-1>; struct { SerializedSCT
 * sct_list <1..2^16-1>; } SignedCertificateTimestampList;
 */
# define MAX_SCT_SIZE            65535
# define MAX_SCT_LIST_SIZE       MAX_SCT_SIZE

/*
 * Macros to read and write integers in network-byte order.
 */

#define n2s(c,s)        ((s=(((unsigned int)((c)[0]))<< 8)| \
                            (((unsigned int)((c)[1]))    )),c+=2)

#define s2n(s,c)        ((c[0]=(unsigned char)(((s)>> 8)&0xff), \
                          c[1]=(unsigned char)(((s)    )&0xff)),c+=2)

#define l2n3(l,c)       ((c[0]=(unsigned char)(((l)>>16)&0xff), \
                          c[1]=(unsigned char)(((l)>> 8)&0xff), \
                          c[2]=(unsigned char)(((l)    )&0xff)),c+=3)

#define n2l8(c,l)       (l =((uint64_t)(*((c)++)))<<56, \
                         l|=((uint64_t)(*((c)++)))<<48, \
                         l|=((uint64_t)(*((c)++)))<<40, \
                         l|=((uint64_t)(*((c)++)))<<32, \
                         l|=((uint64_t)(*((c)++)))<<24, \
                         l|=((uint64_t)(*((c)++)))<<16, \
                         l|=((uint64_t)(*((c)++)))<< 8, \
                         l|=((uint64_t)(*((c)++))))

#define l2n8(l,c)       (*((c)++)=(unsigned char)(((l)>>56)&0xff), \
                         *((c)++)=(unsigned char)(((l)>>48)&0xff), \
                         *((c)++)=(unsigned char)(((l)>>40)&0xff), \
                         *((c)++)=(unsigned char)(((l)>>32)&0xff), \
                         *((c)++)=(unsigned char)(((l)>>24)&0xff), \
                         *((c)++)=(unsigned char)(((l)>>16)&0xff), \
                         *((c)++)=(unsigned char)(((l)>> 8)&0xff), \
                         *((c)++)=(unsigned char)(((l)    )&0xff))

/* Signed Certificate Timestamp */
struct sct_st {
    sct_version_t version;
    /* If version is not SCT_VERSION_V1, this contains the encoded SCT */
    unsigned char *sct;
    size_t sct_len;
    /* If version is SCT_VERSION_V1, fields below contain components of the SCT */
    unsigned char *log_id;
    size_t log_id_len;
    /*
    * Note, we cannot distinguish between an unset timestamp, and one
    * that is set to 0.  However since CT didn't exist in 1970, no real
    * SCT should ever be set as such.
    */
    uint64_t timestamp;
    unsigned char *ext;
    size_t ext_len;
    unsigned char hash_alg;
    unsigned char sig_alg;
    unsigned char *sig;
    size_t sig_len;
    /* Log entry type */
    ct_log_entry_type_t entry_type;
    /* Where this SCT was found, e.g. certificate, OCSP response, etc. */
    sct_source_t source;
    /* The CT log that produced this SCT. */
    CTLOG *log;
    /* The result of the last attempt to validate this SCT. */
    sct_validation_status_t validation_status;
};

/* Miscellaneous data that is useful when verifying an SCT  */
struct sct_ctx_st {
    /* Public key */
    EVP_PKEY *pkey;
    /* Hash of public key */
    unsigned char *pkeyhash;
    size_t pkeyhashlen;
    /* For pre-certificate: issuer public key hash */
    unsigned char *ihash;
    size_t ihashlen;
    /* certificate encoding */
    unsigned char *certder;
    size_t certderlen;
    /* pre-certificate encoding */
    unsigned char *preder;
    size_t prederlen;
};

/* Context when evaluating whether a Certificate Transparency policy is met */
struct ct_policy_eval_ctx_st {
    X509 *cert;
    X509 *issuer;
    CTLOG_STORE *log_store;
    STACK_OF(SCT) *good_scts;
    STACK_OF(SCT) *bad_scts;
};

/*
 * Creates a new context for verifying an SCT.
 */
SCT_CTX *SCT_CTX_new(void);
/*
 * Deletes an SCT verification context.
 */
void SCT_CTX_free(SCT_CTX *sctx);

/* Sets the certificate that the SCT is related to */
int SCT_CTX_set1_cert(SCT_CTX *sctx, X509 *cert, X509 *presigner);
/* Sets the issuer of the certificate that the SCT is related to */
int SCT_CTX_set1_issuer(SCT_CTX *sctx, const X509 *issuer);
/* Sets the public key of the issuer of the certificate that the SCT relates to */
int SCT_CTX_set1_issuer_pubkey(SCT_CTX *sctx, X509_PUBKEY *pubkey);
/* Sets the public key of the CT log that the SCT is from */
int SCT_CTX_set1_pubkey(SCT_CTX *sctx, X509_PUBKEY *pubkey);

/*
 * Does this SCT have the minimum fields populated to be usuable?
 * Returns 1 if so, 0 otherwise.
 */
int SCT_is_complete(const SCT *sct);

/*
 * Does this SCT have the signature-related fields populated?
 * Returns 1 if so, 0 otherwise.
 * This checks that the signature and hash algorithms are set to supported
 * values and that the signature field is set.
 */
int SCT_signature_is_complete(const SCT *sct);


