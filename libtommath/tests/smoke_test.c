#include "tommath.h"

int main(void)
{
   mp_int a, b, c;
   int err;

   if ((err = mp_init_multi(&a, &b, &c, NULL)) != MP_OKAY) {
      return 1;
   }

   mp_set_int(&a, 12345);
   mp_set_int(&b, 67890);
   if ((err = mp_add(&a, &b, &c)) != MP_OKAY) {
      mp_clear_multi(&a, &b, &c, NULL);
      return 1;
   }
   if (mp_cmp_d(&c, 80235) != MP_EQ) {
      mp_clear_multi(&a, &b, &c, NULL);
      return 1;
   }

   mp_clear_multi(&a, &b, &c, NULL);
   return 0;
}
