@with_kw mutable struct TriSolver{FT,N}
  a::Vector{FT} = zeros(FT, N)
  b::Vector{FT} = zeros(FT, N)
  c::Vector{FT} = zeros(FT, N)
  d::Vector{FT} = zeros(FT, N)
  e::Vector{FT} = zeros(FT, N)
  f::Vector{FT} = zeros(FT, N)
  # solution
  u::Vector{FT} = zeros(FT, N)
  du::Vector{FT} = zeros(FT, N)
end

"""
    tridiagonal_solver(a, b, c, d)

Tridiagonal solver

# Description

Converted into a R code from the original code of Gordon Bonan: Bonan, G.
(2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge:
Cambridge University Press. doi:10.1017/9781107339217

Solve for U given the set of equations R * U = D, where U is a vector
of length N, D is a vector of length N, and R is an N x N tridiagonal
matrix defined by the vectors A, B, C each of length N. A(1) and
C(N) are undefined and are not referenced.

    |B(1) C(1) ...  ...  ...                     |
    |A(2) B(2) C(2) ...  ...                     |
R = |     A(3) B(3) C(3) ...                     |
    |                    ... A(N-1) B(N-1) C(N-1)|
    |                    ... ...    A(N)   B(N)  |

The system of equations is written as:

  A_i * U_i-1 + B_i * U_i + C_i * U_i+1 = D_i

for i = 1 to N. The solution is found by rewriting the
equations so that:

  U_i = F_i - E_i * U_i+1

# Return
- `Solution: U`

# Example
```julia
tridiagonal_solver!(a, b, c, d, e, f, u; ibeg)
```
"""
function tridiagonal_solver!(a::V, b::V, c::V, d::V, e::V, f::V, u::V;
  ibeg::Int=1, N::Int=length(a)) where {T<:Real,V<:AbstractVector{T}}

  # e = fill(NaN, N - 1)
  e[ibeg] = c[ibeg] / b[ibeg]
  @inbounds for i in ibeg+1:(N-1)
    e[i] = c[i] / (b[i] - a[i] * e[i-1])
  end

  # f = fill(NaN, N)
  f[ibeg] = d[ibeg] / b[ibeg]
  @inbounds for i in ibeg+1:(N)
    f[i] = (d[i] - a[i] * f[i-1]) / (b[i] - a[i] * e[i-1])
  end

  # u = fill(NaN, N)
  u[N] = f[N]
  @inbounds for i in N-1:-1:ibeg
    u[i] = f[i] - e[i] * u[i+1]
  end
  return u
end

function tridiagonal_solver(a::V, b::V, c::V, d::V;
  ibeg::Int=1, N::Int=length(a)) where {T<:Real,V<:AbstractVector{T}}

  e = fill(NaN, N - 1)
  f = fill(NaN, N)
  u = fill(NaN, N)
  tridiagonal_solver!(a, b, c, d, e, f, u; ibeg, N)
  return u
end


export tridiagonal_solver, tridiagonal_solver!, TriSolver
