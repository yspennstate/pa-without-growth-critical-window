# The critical window of preferential attachment without vertex growth

Code for the paper "The critical window of preferential attachment without vertex growth" (Y. Shmalo).

The paper source and PDF are in `paper/`.

`check.py` uses exact rational arithmetic and does the following:

- It enumerates every path of the simple process and of the multigraph process in four small cases, (n, m, alpha) = (4, 3, 1), (5, 4, 1), (5, 3, 1/4) and (4, 4, 4). On each step it checks the likelihood-ratio identity of Section 2.
- In each of these cases it prints the probability Q_n(A_n) that the multigraph process is simple after m edges and the total variation distance between the simple process and the multigraph process conditioned to be simple. The case n = 5, m = 4, alpha = 1 (5040 paths, Q_n(A_n) = 8/135, distance 122/1463) is the one quoted in the introduction.
- It checks the factorial moments of the negative binomial limit D, and the identity eta = mu^2 kappa.
- It checks that the functions s(t) and r(t) of Lemma 3.2 solve the fluid equations, and that the limit of n b_k is c + 2c^2 t.
- It prints the exact mean of nu - 1 at the critical time for n = 1000.

The output of `python check.py` is in `check_output.txt`. It runs in a few seconds with Python 3 and no other packages.
