/* Quick smoke-test for rat, cusp, and p1 layers. */
#include <stdio.h>
#include <assert.h>
#include "rat.h"
#include "cusp.h"
#include "p1.h"

static void test_rat(void)
{
    Rat a = rat_make(1, 2);
    Rat b = rat_make(1, 3);
    Rat r;

    assert(rat_add(a, b, &r) == RAT_OK);
#if defined(RAT_BACKEND_BUILTIN)
    assert(r.num == 5 && r.den == 6);  /* 1/2 + 1/3 = 5/6 */
#else
    assert(rat_eq(r, rat_make(5, 6)));
#endif

    assert(rat_mul(a, b, &r) == RAT_OK);
#if defined(RAT_BACKEND_BUILTIN)
    assert(r.num == 1 && r.den == 6);  /* 1/2 * 1/3 = 1/6 */
#else
    assert(rat_eq(r, rat_make(1, 6)));
#endif

    assert(rat_sub(a, b, &r) == RAT_OK);
#if defined(RAT_BACKEND_BUILTIN)
    assert(r.num == 1 && r.den == 6);  /* 1/2 - 1/3 = 1/6 */
#else
    assert(rat_eq(r, rat_make(1, 6)));
#endif

    Rat c = rat_make(3, 4);
    Rat d = rat_make(6, 8);
    assert(rat_eq(c, d));              /* both reduce to 3/4 */

    assert(rat_inv(a, &r) == RAT_OK);
#if defined(RAT_BACKEND_BUILTIN)
    assert(r.num == 2 && r.den == 1);  /* inv(1/2) = 2 */
#else
    assert(rat_eq(r, rat_int(2)));
#endif

    /* Division by zero. */
    assert(rat_inv(rat_zero(), &r) == RAT_DIVZERO);

    printf("rat: OK\n");
}

static void test_cusp(void)
{
    Cusp inf = cusp_inf();
    assert(cusp_is_inf(inf));

    Cusp c = cusp_from_frac(4, 6);    /* reduces to 2/3 */
    assert(!cusp_is_inf(c));
    assert(c.p == 2 && c.q == 3);

    Cusp d = cusp_from_frac(-3, 9);   /* reduces to -1/3 */
    assert(d.p == -1 && d.q == 3);

    assert(cusp_eq(cusp_from_frac(1, 2), cusp_from_frac(2, 4)));

    printf("cusp: OK\n");
}

/* Known values: mu = [SL2(Z) : Gamma_0(N)]
 *   N=1:  mu=1
 *   N=2:  mu=3
 *   N=3:  mu=4
 *   N=4:  mu=6
 *   N=5:  mu=6
 *   N=6:  mu=12
 *   N=11: mu=12  */
static void test_p1_mu(void)
{
    int expected[][2] = {{1,1},{2,3},{3,4},{4,6},{5,6},{6,12},{11,12}};
    int n = sizeof(expected) / sizeof(expected[0]);
    for (int i = 0; i < n; i++) {
        int N = expected[i][0], mu = expected[i][1];
        P1Table *t = p1table_build(N);
        assert(t != NULL);
        if (t->mu != mu) {
            printf("FAIL: P1(Z/%dZ) mu=%d, expected %d\n", N, t->mu, mu);
        } else {
            printf("P1(Z/%dZ): mu=%d OK\n", N, t->mu);
        }
        p1table_free(t);
    }
}

static void test_p1_lift(void)
{
    /* Verify that every element lifts to a valid SL2(Z) matrix. */
    int Ns[] = {2, 3, 5, 11};
    for (int ni = 0; ni < 4; ni++) {
        int N = Ns[ni];
        P1Table *t = p1table_build(N);
        for (int i = 0; i < t->mu; i++) {
            int64_t a, b, c, d;
            int rc = p1_lift_to_sl2z(t, i, &a, &b, &c, &d);
            assert(rc == 0);
            assert(a * d - b * c == 1);   /* det == 1 */
        }
        printf("p1_lift N=%d: OK\n", N);
        p1table_free(t);
    }
}

static void test_p1_S_involution(void)
{
    /* S has order 2 in PSL2(Z): applying S twice returns to start. */
    P1Table *t = p1table_build(11);
    for (int i = 0; i < t->mu; i++) {
        int si = p1_apply_S(t, i);
        assert(si >= 0);
        int ssi = p1_apply_S(t, si);
        assert(ssi == i);
    }
    printf("p1_apply_S involution N=11: OK\n");
    p1table_free(t);
}

static void test_p1_U_order3(void)
{
    /* U has order 3 in PSL2(Z): applying U three times returns to start. */
    P1Table *t = p1table_build(11);
    for (int i = 0; i < t->mu; i++) {
        int ui   = p1_apply_U(t, i);
        int uui  = p1_apply_U(t, ui);
        int uuui = p1_apply_U(t, uui);
        assert(ui >= 0 && uui >= 0 && uuui == i);
    }
    printf("p1_apply_U order-3 N=11: OK\n");
    p1table_free(t);
}

int main(void)
{
    test_rat();
    test_cusp();
    test_p1_mu();
    test_p1_lift();
    test_p1_S_involution();
    test_p1_U_order3();
    printf("\nAll tests passed.\n");
    return 0;
}
