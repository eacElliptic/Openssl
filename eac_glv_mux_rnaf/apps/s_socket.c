/* Copyright (C) 1995-1998 Eric Young (eay@cryptsoft.com)
 * All rights reserved.
 *
 * This package is an SSL implementation written
 * by Eric Young (eay@cryptsoft.com).
 * The implementation was written so as to conform with Netscapes SSL.
 *
 * This library is free for commercial and non-commercial use as long as
 * the following conditions are aheared to.  The following conditions
 * apply to all code found in this distribution, be it the RC4, RSA,
 * lhash, DES, etc., code; not just the SSL code.  The SSL documentation
 * included with this distribution is covered by the same copyright terms
 * except that the holder is Tim Hudson (tjh@cryptsoft.com).
 *
 * Copyright remains Eric Young's, and as such any Copyright notices in
 * the code are not to be removed.
 * If this package is used in a product, Eric Young should be given attribution
 * as the author of the parts of the library used.
 * This can be in the form of a textual message at program startup or
 * in documentation (online or textual) provided with the package.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    "This product includes cryptographic software written by
 *     Eric Young (eay@cryptsoft.com)"
 *    The word 'cryptographic' can be left out if the rouines from the library
 *    being used are not cryptographic related :-).
 * 4. If you include any Windows specific code (or a derivative thereof) from
 *    the apps directory (application code) you must include an acknowledgement:
 *    "This product includes software written by Tim Hudson (tjh@cryptsoft.com)"
 *
 * THIS SOFTWARE IS PROVIDED BY ERIC YOUNG ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * The licence and distribution terms for any publically available version or
 * derivative of this code cannot be changed.  i.e. this code cannot simply be
 * copied and put under another distribution licence
 * [including the GNU Public Licence.]
 */
/* ====================================================================
 * Copyright (c) 199-2015 The OpenSSL Project.  All rights reserved.
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

/* socket-related functions used by s_client and s_server */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h>

/*
 * With IPv6, it looks like Digital has mixed up the proper order of
 * recursive header file inclusion, resulting in the compiler complaining
 * that u_int isn't defined, but only if _POSIX_C_SOURCE is defined, which is
 * needed to have fileno() declared correctly...  So let's define u_int
 */
#if defined(OPENSSL_SYS_VMS_DECC) && !defined(__U_INT)
# define __U_INT
typedef unsigned int u_int;
#endif

#ifndef OPENSSL_NO_SOCK

# define USE_SOCKETS
# include "apps.h"
# undef USE_SOCKETS
# include "s_apps.h"

# include <openssl/bio.h>
# include <openssl/err.h>

/*
 * init_client - helper routine to set up socket communication
 * @sock: pointer to storage of resulting socket.
 * @host: the host name or path (for AF_UNIX) to connect to.
 * @port: the port to connect to (ignored for AF_UNIX).
 * @family: desired socket family, may be AF_INET, AF_INET6, AF_UNIX or
 *  AF_UNSPEC
 * @type: socket type, must be SOCK_STREAM or SOCK_DGRAM
 *
 * This will create a socket and use it to connect to a host:port, or if
 * family == AF_UNIX, to the path found in host.
 *
 * If the host has more than one address, it will try them one by one until
 * a successful connection is established.  The resulting socket will be
 * found in *sock on success, it will be given INVALID_SOCKET otherwise.
 *
 * Returns 1 on success, 0 on failure.
 */
int init_client(int *sock, const char *host, const char *port,
                int family, int type)
{
    BIO_ADDRINFO *res = NULL;
    const BIO_ADDRINFO *ai = NULL;
    int ret;

    if (!BIO_sock_init())
        return 0;

    ret = BIO_lookup(host, port, BIO_LOOKUP_CLIENT, family, type, &res);
    if (ret == 0) {
        ERR_print_errors(bio_err);
        return 0;
    }

    ret = 0;
    for (ai = res; ai != NULL; ai = BIO_ADDRINFO_next(ai)) {
        /* Admittedly, these checks are quite paranoid, we should not get
         * anything in the BIO_ADDRINFO chain that we haven't
         * asked for. */
        OPENSSL_assert((family == AF_UNSPEC || family == BIO_ADDRINFO_family(res))
                       && (type == 0 || type == BIO_ADDRINFO_socktype(res)));

        *sock = BIO_socket(BIO_ADDRINFO_family(ai), BIO_ADDRINFO_socktype(ai),
                           BIO_ADDRINFO_protocol(res), 0);
        if (*sock == INVALID_SOCKET) {
            /* Maybe the kernel doesn't support the socket family, even if
             * BIO_lookup() added it in the returned result...
             */
            continue;
        }
        if (!BIO_connect(*sock, BIO_ADDRINFO_address(ai), 0)) {
            BIO_closesocket(*sock);
            *sock = INVALID_SOCKET;
            continue;
        }

        /* Success, don't try any more addresses */
        break;
    }

    if (*sock == INVALID_SOCKET) {
        ERR_print_errors(bio_err);
    } else {
        ret = 1;
    }
    BIO_ADDRINFO_free(res);
    return ret;
}

/*
 * do_server - helper routine to perform a server operation
 * @accept_sock: pointer to storage of resulting socket.
 * @host: the host name or path (for AF_UNIX) to connect to.
 * @port: the port to connect to (ignored for AF_UNIX).
 * @family: desired socket family, may be AF_INET, AF_INET6, AF_UNIX or
 *  AF_UNSPEC
 * @type: socket type, must be SOCK_STREAM or SOCK_DGRAM
 * @cb: pointer to a function that receives the accepted socket and
 *  should perform the communication with the connecting client.
 * @context: pointer to memory that's passed verbatim to the cb function.
 * @naccept: number of times an incoming connect should be accepted.  If -1,
 *  unlimited number.
 *
 * This will create a socket and use it to listen to a host:port, or if
 * family == AF_UNIX, to the path found in host, then start accepting
 * incoming connections and run cb on the resulting socket.
 *
 * 0 on failure, something other on success.
 */
int do_server(int *accept_sock, const char *host, const char *port,
              int family, int type, do_server_cb cb,
              unsigned char *context, int naccept)
{
    int asock = 0;
    int sock;
    int i;
    BIO_ADDRINFO *res = NULL;
    int ret = 0;

    if (!BIO_sock_init())
        return 0;

    if (!BIO_lookup(host, port, BIO_LOOKUP_SERVER, family, type, &res)) {
        ERR_print_errors(bio_err);
        return 0;
    }

    /* Admittedly, these checks are quite paranoid, we should not get
     * anything in the BIO_ADDRINFO chain that we haven't asked for */
    OPENSSL_assert((family == AF_UNSPEC || family == BIO_ADDRINFO_family(res))
                   && (type == 0 || type == BIO_ADDRINFO_socktype(res)));

    asock = BIO_socket(BIO_ADDRINFO_family(res), BIO_ADDRINFO_socktype(res),
                       BIO_ADDRINFO_protocol(res), 0);
    if (asock == INVALID_SOCKET
        || !BIO_listen(asock, BIO_ADDRINFO_address(res), BIO_SOCK_REUSEADDR)) {
        BIO_ADDRINFO_free(res);
        ERR_print_errors(bio_err);
        if (asock != INVALID_SOCKET)
            BIO_closesocket(asock);
        goto end;
    }

    BIO_ADDRINFO_free(res);
    res = NULL;

    if (accept_sock != NULL)
        *accept_sock = asock;
    for (;;) {
        if (type == SOCK_STREAM) {
            do {
                sock = BIO_accept_ex(asock, NULL, 0);
            } while (sock < 0 && BIO_sock_should_retry(ret));
            if (sock < 0) {
                ERR_print_errors(bio_err);
                SHUTDOWN(asock);
                break;
            }
            i = (*cb)(sock, type, context);
            SHUTDOWN2(sock);
        } else {
            i = (*cb)(asock, type, context);
        }

        if (naccept != -1)
            naccept--;
        if (i < 0 || naccept == 0) {
            SHUTDOWN2(asock);
            ret = i;
            break;
        }
    }
 end:
# ifdef AF_UNIX
    if (family == AF_UNIX)
        unlink(host);
# endif
    return ret;
}

#endif  /* OPENSSL_NO_SOCK */