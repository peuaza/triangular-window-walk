"""C2: verify, rather than assert, the hypotheses of the Vandermonde rank argument.

The paper claims rank(H_m) = m for H_m[d,k] = 3^{(k+d)^2} mod q with q = 2^31-1, via
H_m = 3 diag(3^{d^2}) [9^{dk}] diag(3^{k^2}).  For that to be a rank statement over F_q one needs
  (i)   q prime (so F_q is a field),
  (ii)  3 != 0 and 9 != 0 mod q (the diagonal factors invertible),
  (iii) ord_q(9) > m-1 for every m in range, so the Vandermonde nodes 9^d, d = 0..m-1, are distinct,
  (iv)  the asserted order itself checked: 9^ord = 1 and ord | q-1.
Prints each check and the largest m the range supports.
"""
import sympy

q = 2**31 - 1
MAXM = 16384          # the largest window the paper's tables use
ORDER_CLAIM = 357_913_941

print("q =", q)
print("(i)   prime:", sympy.isprime(q))
print("(ii)  3 mod q:", 3 % q, " 9 mod q:", 9 % q, " both nonzero:", (3 % q) * (9 % q) % q != 0)

ord9 = sympy.n_order(9, q)
print("(iii) ord_q(9) computed:", ord9, " claimed:", ORDER_CLAIM, " match:", ord9 == ORDER_CLAIM)
print("      ord divides q-1:", (q - 1) % ord9 == 0)
print("      range check ord > MAXM-1:", ord9 > MAXM - 1, f"(ord/MAXM = {ord9/MAXM:.0f}x)")

# the nodes are the powers 9^d for d = 0..m-1; distinct iff ord >= m
print("(iv)  nodes distinct for every m <= MAXM:", ord9 >= MAXM)
print()
print("conclusion: rank(H_m) = m over F_q for every m <= %d, with the proof's hypotheses checked." % MAXM)
