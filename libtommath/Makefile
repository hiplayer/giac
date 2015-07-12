include ../build/core.mk
LOCAL_PATH:= $(call my-dir)
include $(CLEAR_VARS)

LOCAL_MODULE := libtommath

LOCAL_SRC_FILES:= \
	libtommath-0.39/bncore.c \
	libtommath-0.39/bn_mp_init.c \
	libtommath-0.39/bn_mp_clear.c \
	libtommath-0.39/bn_mp_exch.c \
	libtommath-0.39/bn_mp_grow.c \
	libtommath-0.39/bn_mp_shrink.c \
	libtommath-0.39/bn_mp_clamp.c \
	libtommath-0.39/bn_mp_zero.c \
	libtommath-0.39/bn_mp_set.c \
	libtommath-0.39/bn_mp_set_int.c \
	libtommath-0.39/bn_mp_init_size.c \
	libtommath-0.39/bn_mp_copy.c \
	libtommath-0.39/bn_mp_init_copy.c \
	libtommath-0.39/bn_mp_abs.c \
	libtommath-0.39/bn_mp_neg.c \
	libtommath-0.39/bn_mp_cmp_mag.c \
	libtommath-0.39/bn_mp_cmp.c \
	libtommath-0.39/bn_mp_cmp_d.c \
	libtommath-0.39/bn_mp_rshd.c \
	libtommath-0.39/bn_mp_lshd.c \
	libtommath-0.39/bn_mp_mod_2d.c \
	libtommath-0.39/bn_mp_div_2d.c \
	libtommath-0.39/bn_mp_mul_2d.c \
	libtommath-0.39/bn_mp_div_2.c \
	libtommath-0.39/bn_mp_mul_2.c \
	libtommath-0.39/bn_s_mp_add.c \
	libtommath-0.39/bn_s_mp_sub.c \
	libtommath-0.39/bn_fast_s_mp_mul_digs.c \
	libtommath-0.39/bn_s_mp_mul_digs.c \
	libtommath-0.39/bn_fast_s_mp_mul_high_digs.c \
	libtommath-0.39/bn_s_mp_mul_high_digs.c \
	libtommath-0.39/bn_fast_s_mp_sqr.c \
	libtommath-0.39/bn_s_mp_sqr.c \
	libtommath-0.39/bn_mp_add.c \
	libtommath-0.39/bn_mp_sub.c \
	libtommath-0.39/bn_mp_karatsuba_mul.c \
	libtommath-0.39/bn_mp_mul.c \
	libtommath-0.39/bn_mp_karatsuba_sqr.c \
	libtommath-0.39/bn_mp_sqr.c \
	libtommath-0.39/bn_mp_div.c \
	libtommath-0.39/bn_mp_mod.c \
	libtommath-0.39/bn_mp_add_d.c \
	libtommath-0.39/bn_mp_sub_d.c \
	libtommath-0.39/bn_mp_mul_d.c \
	libtommath-0.39/bn_mp_div_d.c \
	libtommath-0.39/bn_mp_mod_d.c \
	libtommath-0.39/bn_mp_expt_d.c \
	libtommath-0.39/bn_mp_addmod.c \
	libtommath-0.39/bn_mp_submod.c \
	libtommath-0.39/bn_mp_mulmod.c \
	libtommath-0.39/bn_mp_sqrmod.c \
	libtommath-0.39/bn_mp_gcd.c \
	libtommath-0.39/bn_mp_lcm.c \
	libtommath-0.39/bn_fast_mp_invmod.c \
	libtommath-0.39/bn_mp_invmod.c \
	libtommath-0.39/bn_mp_reduce.c \
	libtommath-0.39/bn_mp_montgomery_setup.c \
	libtommath-0.39/bn_fast_mp_montgomery_reduce.c \
	libtommath-0.39/bn_mp_montgomery_reduce.c \
	libtommath-0.39/bn_mp_exptmod_fast.c \
	libtommath-0.39/bn_mp_exptmod.c \
	libtommath-0.39/bn_mp_2expt.c \
	libtommath-0.39/bn_mp_n_root.c \
	libtommath-0.39/bn_mp_jacobi.c \
	libtommath-0.39/bn_reverse.c \
	libtommath-0.39/bn_mp_count_bits.c \
	libtommath-0.39/bn_mp_read_unsigned_bin.c \
	libtommath-0.39/bn_mp_read_signed_bin.c \
	libtommath-0.39/bn_mp_to_unsigned_bin.c \
	libtommath-0.39/bn_mp_to_signed_bin.c \
	libtommath-0.39/bn_mp_unsigned_bin_size.c \
	libtommath-0.39/bn_mp_signed_bin_size.c \
	libtommath-0.39/bn_mp_xor.c \
	libtommath-0.39/bn_mp_and.c \
	libtommath-0.39/bn_mp_or.c \
	libtommath-0.39/bn_mp_rand.c \
	libtommath-0.39/bn_mp_montgomery_calc_normalization.c \
	libtommath-0.39/bn_mp_prime_is_divisible.c \
	libtommath-0.39/bn_prime_tab.c \
	libtommath-0.39/bn_mp_prime_fermat.c \
	libtommath-0.39/bn_mp_prime_miller_rabin.c \
	libtommath-0.39/bn_mp_prime_is_prime.c \
	libtommath-0.39/bn_mp_prime_next_prime.c \
	libtommath-0.39/bn_mp_dr_reduce.c \
	libtommath-0.39/bn_mp_dr_is_modulus.c \
	libtommath-0.39/bn_mp_dr_setup.c \
	libtommath-0.39/bn_mp_reduce_setup.c \
	libtommath-0.39/bn_mp_toom_mul.c \
	libtommath-0.39/bn_mp_toom_sqr.c \
	libtommath-0.39/bn_mp_div_3.c \
	libtommath-0.39/bn_s_mp_exptmod.c \
	libtommath-0.39/bn_mp_reduce_2k.c \
	libtommath-0.39/bn_mp_reduce_is_2k.c \
	libtommath-0.39/bn_mp_reduce_2k_setup.c \
	libtommath-0.39/bn_mp_reduce_2k_l.c \
	libtommath-0.39/bn_mp_reduce_is_2k_l.c \
	libtommath-0.39/bn_mp_reduce_2k_setup_l.c \
	libtommath-0.39/bn_mp_radix_smap.c \
	libtommath-0.39/bn_mp_read_radix.c \
	libtommath-0.39/bn_mp_toradix.c \
	libtommath-0.39/bn_mp_radix_size.c \
	libtommath-0.39/bn_mp_fread.c \
	libtommath-0.39/bn_mp_fwrite.c \
	libtommath-0.39/bn_mp_cnt_lsb.c \
	libtommath-0.39/bn_error.c \
	libtommath-0.39/bn_mp_init_multi.c \
	libtommath-0.39/bn_mp_clear_multi.c \
	libtommath-0.39/bn_mp_exteuclid.c \
	libtommath-0.39/bn_mp_toradix_n.c \
	libtommath-0.39/bn_mp_prime_random_ex.c \
	libtommath-0.39/bn_mp_get_int.c \
	libtommath-0.39/bn_mp_sqrt.c \
	libtommath-0.39/bn_mp_is_square.c \
	libtommath-0.39/bn_mp_init_set.c \
	libtommath-0.39/bn_mp_init_set_int.c \
	libtommath-0.39/bn_mp_invmod_slow.c \
	libtommath-0.39/bn_mp_prime_rabin_miller_trials.c \
	libtommath-0.39/bn_mp_to_signed_bin_n.c \
	libtommath-0.39/bn_mp_to_unsigned_bin_n.c \

LOCAL_CFLAGS := \
	-Wall -W -Wshadow -Wsign-compare \

LOCAL_C_INCLUDES := \
	$(LOCAL_CFLAGS) \
	-Ilibtommath-0.39 \
	-I../include \

LOCAL_SHARED_LIBRARIES := \

LOCAL_RELEASE_TYPE := debug

include $(BUILD_STATIC_LIBRARY)
