#include <openssl/err.h>
#include "ec_lcl.h"

/**
 * Copyright OpenSSL 2016
 * Contents licensed under the terms of the OpenSSL license
 * See http://www.openssl.org/source/license.html for details
 *
 * Euclidean addition chains scalar multiplication on curves with efficient endomorphism
 */



/* GLV-related per-curve constants */
 static const unsigned char eac_beta_c1p358k256[] = {
    /* beta */
    0x34, 0x24, 0x4B, 0xA4, 0xDB, 0xD4, 0xE3, 0x52, 0x73, 0xBE, 0xBD, 0x1E, 0x4E, 0x1D,
    0x25, 0xA1, 0x8A, 0x5A, 0x79, 0xAB, 0xAE, 0x77, 0xF4, 0x3C, 0xE1, 0x88, 0x89, 0x68,
    0x2A, 0xBB, 0xF5, 0xFE, 0x04, 0xBE, 0xB4, 0xAB, 0xAD, 0x9E, 0xE0, 0xDA, 0xB9, 0x28,
    0x98, 0xEA, 0x54
};



typedef struct eac_data{
	int eac_length ;
	BIGNUM *beta;
} eac_data;




int ec_GFp_eac_group_init(EC_GROUP *group){
    int ok;

    ok = ec_GFp_mont_group_init(group);
    group->custom_data = NULL;
    return ok;
}


#define EAC_CONSTANTS_FREE(g) do {                   \
    if (g->custom_data != NULL) {                    \
        eac_data *data = (eac_data *)g->custom_data; \
        BN_free(data->beta);                   		 \                
        OPENSSL_free(g->custom_data);                \
        g->custom_data = NULL;                       \
    }                                                \
} while(0)


void copy_data(eac_data *dest, eac_data *src){
	dest->eac_length = src->eac_length;
	dest->beta = BN_dup(src->beta);
}


void ec_GFp_eac_group_finish(EC_GROUP *group)
{
    EAC_CONSTANTS_FREE(group);
    ec_GFp_mont_group_finish(group);
}


void ec_GFp_eac_group_clear_finish(EC_GROUP *group)
{

    if (group->custom_data != NULL) {
        eac_data *data = (eac_data *)group->custom_data;
        BN_clear_free(data->beta);
        OPENSSL_clear_free(group->custom_data, sizeof(eac_data));
        group->custom_data = NULL;
    }
    ec_GFp_mont_group_clear_finish(group);
}


int ec_GFp_eac_group_copy(EC_GROUP *dest, const EC_GROUP *src)
{

    if (!ec_GFp_mont_group_copy(dest, src))
        return 0;

    EAC_CONSTANTS_FREE(dest);

    if (src->custom_data == NULL)
        return 0;

    dest->custom_data = OPENSSL_zalloc(sizeof(eac_data));
    if (dest->custom_data == NULL)
        return 0;
    
    eac_data *destc = (eac_data *)dest->custom_data;
    eac_data *srcc = (eac_data *)src->custom_data;
    
    copy_data(destc, srcc);
    
    return 1;
}


int ec_GFp_eac_group_set_curve(EC_GROUP *group, const BIGNUM *p,
                               const BIGNUM *a, const BIGNUM *b, BN_CTX *ctx)
{
    if (ec_GFp_mont_group_set_curve(group, p, a, b, ctx) == 0)
        return 0;

    EAC_CONSTANTS_FREE(group);

    group->custom_data = OPENSSL_zalloc(sizeof(eac_data));
    if (group->custom_data == NULL)
        return 0;
    
    eac_data *data = (eac_data *)group->custom_data;
    
	if ((data->beta = BN_new()) == NULL)
		goto err;
		
    switch (group->curve_name) {
		case NID_eac256k:
			data->eac_length = 256;
			BN_bin2bn(eac_beta_c1p358k256, 45, data->beta);
			break;
		default:
			goto err;
    }

	if (data->beta == NULL)
		goto err;

	
	/* encode beta parameter to curve's finite field */
    if (!group->meth->field_encode(group, data->beta, data->beta, ctx))
        goto err;
	
    return 1;
    

 err:
    EAC_CONSTANTS_FREE(group);
    return 0;
}



//~ Always computes (p1, p1+p2) with p1 updated.
int ZADDU(const EC_GROUP *group, EC_POINT *p1, EC_POINT *p2, BIGNUM *A, BIGNUM *B, BN_CTX *ctx){
	//~ A <- (X2-X1)^2
	if (!group->meth->field_sqr(group, A, B, ctx))
		return 0;
	//~ X1 <- X1*(X2-X1)^2
	if (!group->meth->field_mul(group, p1->X, p1->X, A, ctx))
		return 0;
	//~ A <- X2*(X2-X1)^2 
	if (!group->meth->field_mul(group, A, p2->X, A, ctx))
		return 0;
	//~ B <- (Y2-Y1)^2 
	if (!BN_sub(p2->Y, p2->Y, p1->Y))							// Y2 <- (Y2-Y1)
		return 0;							
	if (!group->meth->field_sqr(group, B, p2->Y, ctx))			// (Y2-Y1)^2 	
		return 0;			
	//~ X2 <- B - X1 - A
	if (!BN_sub(B, B, p1->X))
		return 0;
	if (!BN_mod_sub(p2->X, B, A, group->field, ctx))
		return 0;
	//~ Y1 <- Y1*(A - X1) = Y1*(X2-X1)^3
	if (!BN_sub(A, A, p1->X))									// A <- (A-X1)
		return 0;									
	if (!group->meth->field_mul(group, p1->Y, p1->Y, A, ctx))
		return 0;
	//~ Y2 <- Y2*(X1-X2)-Y1
	if (!BN_sub(B, p1->X, p2->X))								// B <- (X1-X2)
		return 0;								
	if (!group->meth->field_mul(group, B, B, p2->Y, ctx))
		return 0;
	if (!BN_mod_sub_quick(p2->Y, B, p1->Y, group->field))
		return 0;
		
	return 1;
}


int point_from_eac(const EC_GROUP *group, EC_POINT **p0p1, const BIGNUM *eac_k, int eac_len, BIGNUM *A, BIGNUM *B, BN_CTX *ctx){
	int j, n;
		
	//~ Here, we consider (and update) only Z-cordonate of p0p1[1]
	for(j=0; j<eac_len; j++){
		n = BN_is_bit_set(eac_k, j);
		
		//~ B <- X2-X1 
		if (!BN_sub(B, (*(p0p1 + n))->X, (*(p0p1 + (1-n)))->X))
			return 0;
		//~ Z <- Z(X2-X1) 
		if (!group->meth->field_mul(group, p0p1[1]->Z, p0p1[1]->Z, B, ctx))
			return 0;	
		
		if (!ZADDU(group, *(p0p1 + (1-n)), *(p0p1 + n), A, B, ctx))
			return 0;
	}
	
	//~ B <- X2-X1 
	if (!BN_sub(B, (*(p0p1 + 1))->X, (*(p0p1))->X))
		return 0;
	//~ Z <- Z(X2-X1)
	if (!group->meth->field_mul(group, p0p1[1]->Z, p0p1[1]->Z, B, ctx))
		return 0;
	
	p0p1[1]->Z_is_one = BN_is_one(p0p1[1]->Z);    /* useless in the loop 'for' above */
	
	return ZADDU(group, *p0p1, *(p0p1+1), A, B, ctx);
}


/**
 * Computes the sum
 * scalar*group->generator + scalars[0]*points[0] + ... + scalars[num-1]*points[num-1]
 *
 * Here, scalar(s) represent eac.
 */
int ec_GFp_eac_mul(const EC_GROUP *group, EC_POINT *r, const BIGNUM *scalar, size_t num, const EC_POINT *points[], const BIGNUM *scalars[], BN_CTX *ctx){
	
	int i, ret=0;
    eac_data *data;
    EC_POINT **tpoints = NULL;
    BIGNUM *A;
	BIGNUM *B;
	BIGNUM *beta;

	data = (eac_data *)group->custom_data;
    beta = data->beta;

    BN_CTX_start(ctx);
    
    A = BN_CTX_get(ctx);
    B = BN_CTX_get(ctx);
    
    
	/* setup some arrays, apply an endomorphism if scalar is present */
	if (scalar == NULL) {
        if ((tpoints = OPENSSL_malloc((2 * num) * sizeof(EC_POINT *))) == NULL)
            goto end;
    } else {
        if ((tpoints = OPENSSL_malloc((2 * num + 2) * sizeof(EC_POINT *))) == NULL)
            goto end;
        if ((tpoints[2 * num] = EC_POINT_dup(EC_GROUP_get0_generator(group), group)) == NULL)
            goto end;
		if ((tpoints[2 * num + 1] = EC_POINT_dup(EC_GROUP_get0_generator(group), group)) == NULL)
            goto end;
        if (!group->meth->field_mul(group, tpoints[2 * num + 1]->X, tpoints[2 * num + 1]->X, beta, ctx))
            goto end;
    }

	/* apply the endomorphism on all the other points */
    for (i = 0; i < num; i++) {
        if((tpoints[2 * i] = EC_POINT_dup(*((EC_POINT **)points + i), group)) == NULL)
			goto end;
        if((tpoints[2 * i + 1] = EC_POINT_dup(*((EC_POINT **)points + i), group)) == NULL)
            goto end;
        if(!group->meth->field_mul(group, tpoints[2 * i + 1]->X, tpoints[2 * i + 1]->X, beta, ctx))
            goto end;
    }
	
	
	/* result setup */
	 if(!EC_POINT_set_to_infinity(group, r))
		goto end;
	
	/* multiplications */
	if (scalar != NULL) {
		if (!point_from_eac(group, (tpoints + 2 * num), scalar, data->eac_length, A, B, ctx))
			goto end;
		if (!EC_POINT_add(group, r, r, *(tpoints + 2 * num + 1), ctx))
			goto end;
    }
	for (i = 0; i < num; i++) {
		if (!point_from_eac(group, (tpoints + 2 * i), scalars[i], data->eac_length, A, B, ctx))
			goto end;
		if (!EC_POINT_add(group, r, r, *(tpoints + 2 * i + 1), ctx))
			goto end;
	}


	/* everything is ok */
    ret = 1;
    

 end:
	if (tpoints != NULL) {
		for (i = 0; i < num; i++) {
			EC_POINT_free(tpoints[2 * i]);
			EC_POINT_free(tpoints[2 * i + 1]);
		}
		if (scalar != NULL) {
			EC_POINT_free(tpoints[2 * num]);
			EC_POINT_free(tpoints[2 * num + 1]);
		}
	}
	OPENSSL_free(tpoints);
	
	BN_CTX_end(ctx);
	
	return ret;
}




int ec_key_eac_generate_key(EC_KEY *eckey){
	
    int ok = 0, privK_length;
    BN_CTX *ctx = NULL;
    BIGNUM *priv_key = NULL;
    EC_POINT *pub_key = NULL;

    if ((ctx = BN_CTX_new()) == NULL)
        goto err;

    if (eckey->priv_key == NULL) {
        priv_key = BN_new();
        if (priv_key == NULL)
            goto err;
    } else
        priv_key = eckey->priv_key;

    
    privK_length = ((eac_data *)(eckey->group->custom_data))->eac_length;
	if (!BN_rand(priv_key, privK_length, -1, 0))  /* (-1 -> BN_RAND_TOP_ANY) ,  (0 -> BN_RAND_BOTTOM_ANY) */
		goto err;


    if (eckey->pub_key == NULL) {
        pub_key = EC_POINT_new(eckey->group);
        if (pub_key == NULL)
            goto err;
    } else
        pub_key = eckey->pub_key;

    if (!EC_POINT_mul(eckey->group, pub_key, priv_key, NULL, NULL, ctx))
        goto err;


    eckey->priv_key = priv_key;
    eckey->pub_key = pub_key;
    

    ok = 1;


 err:
    if (eckey->pub_key == NULL)
        EC_POINT_free(pub_key);
    if (eckey->priv_key != priv_key)
        BN_free(priv_key);
    BN_CTX_free(ctx);
    return ok;
}




const EC_METHOD *EC_GFp_eac_method(void){
	
    static const EC_METHOD ret = {
        EC_FLAGS_DEFAULT_OCT,
        NID_X9_62_prime_field,
        ec_GFp_eac_group_init,									// new
        ec_GFp_eac_group_finish,								// new
        ec_GFp_eac_group_clear_finish,							// new
        ec_GFp_eac_group_copy,									// new
        ec_GFp_eac_group_set_curve,								// new
        ec_GFp_simple_group_get_curve,
        ec_GFp_simple_group_get_degree,
        ec_group_simple_order_bits,
        ec_GFp_simple_group_check_discriminant,
        ec_GFp_simple_point_init,
        ec_GFp_simple_point_finish,
        ec_GFp_simple_point_clear_finish,
        ec_GFp_simple_point_copy,
        ec_GFp_simple_point_set_to_infinity,
        ec_GFp_simple_set_Jprojective_coordinates_GFp,
        ec_GFp_simple_get_Jprojective_coordinates_GFp,
        ec_GFp_simple_point_set_affine_coordinates,
        ec_GFp_simple_point_get_affine_coordinates,
        0, 0, 0,
        ec_GFp_simple_add,
        ec_GFp_simple_dbl,
        ec_GFp_simple_invert,
        ec_GFp_simple_is_at_infinity,
        ec_GFp_simple_is_on_curve,
        ec_GFp_simple_cmp,
        ec_GFp_simple_make_affine,
        ec_GFp_simple_points_make_affine,
        ec_GFp_eac_mul,											// new
        0 /* precompute_mult */ ,
        0 /* have_precompute_mult */ ,
        ec_GFp_mont_field_mul,
        ec_GFp_mont_field_sqr,
        0, /* field_div */ 
        ec_GFp_mont_field_encode, 
        ec_GFp_mont_field_decode, 
        ec_GFp_mont_field_set_to_one, 
        ec_key_simple_priv2oct,
        ec_key_simple_oct2priv,
        0, /* set private */
        ec_key_eac_generate_key,								// new
        ec_key_simple_check_key,   
        ec_key_simple_generate_public_key, 
        0, /* keycopy */
        0, /* keyfinish */
        ecdh_simple_compute_key                               
    };

    return &ret;
}

