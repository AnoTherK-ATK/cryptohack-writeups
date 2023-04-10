## Quadratic Residues
### Solution
You can use brute force to solve this.

### Code

```python
p = 29
ints = [14, 6, 11]

ans = [x for x in range(p) if(pow(x, 2, p) in ints)]
print(min(ans))
```

## Lengendre Symbol
### Solution
First, we have to find the quadratic residue
`(a/p) = a^((p + 1)/2) mod p`
If `(a/p) = 1` then `a` is quadratic residue
Because `p = 3 mod 4`  so to calculate square root of `a`, we can use [Tonelli–Shanks_algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm)
`a = x^((p + 1) / 4) mod p`

### Code

```python
p = #please give p manually because it's too long

ints = #please give ints manually because it's too long

quad = [x for x in ints if(pow(x, (p - 1)//2, p) == 1)]
# there is one quadratic residue
res = quad[0]
print(pow(res, (p + 1)//4, p))
```

## Modular Square Root
### Solution
Use [Tonelli–Shanks_algorithm](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm) to solve

### Code

```python
def legendre(a, p):
    return pow(a, (p - 1) // 2, p)

def tonelli(n, p):
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r

a = #please give a manually
p = #please give p manually

root = tonelli(a, p)
print(min(root, p - root))
```

## Chinese Remainder Theorem
### Solution
We can calculate with [this algorithm](https://www.geeksforgeeks.org/introduction-to-chinese-remainder-theorem/)

### Code

```python
from Crypto.Util.number import inverse

def crt(a, m):
    Mul = 1
    for i in m:
        Mul *= i
    M = [Mul // x for x in m]
    y = [inverse(M[i], m[i]) for i in range(len(m))]
    ans = 0
    for i in range(len(m)):
        ans += a[i] * M[i] * y[i]
    return ans % Mul

a = [2, 3, 5]
m = [5, 11, 17]
print(crt(a, m))
```

## Vectors
### Solution
You can use `numpy.array()` to solve this

### Code

```python
import numpy as np
v = np.array([2, 6, 3])
w = np.array([1, 0, 0])
u = np.array([7, 7, 2])

x = 3*(2*v - w)
y = 2*u

print(x.dot(y))
```

## Size and Basis
### Solution
Yes, `numpy.array()` too.

### Code

```python
import numpy as np
v = np.array([4, 6, 2, 5])

print(pow(v.dot(v), 0.5))
```

## Gram Schmidt
### Solution
Just use `numpy.array()` and Gram-Schimidt algorithm presented in the challenge

### Code

```python
import numpy as np
v = [np.array([4,1,3,-1]),
     np.array([2,1,-3,4]),
     np.array([1,0,-2,7]),
     np.array([6,2,9,-5])]

u = [v[0]]
for i in range(1, 4):
    mi = [np.dot(v[i], u[j]) / np.dot(u[j], u[j]) for j in range(len(u))]
    u += [v[i] - sum([mij * uj for (mij, uj) in zip(mi, u)])]

print(round(u[3][1], 5))
```

## What's a Lattice
### Solution
Use `numpy.linalg.det()` to solve

### Code

```python
import numpy as np
v = np.array([[6, 2, -3], [5, 1, 4], [2, 7, 1]])

print(round(abs(np.linalg.det(v))))
```
