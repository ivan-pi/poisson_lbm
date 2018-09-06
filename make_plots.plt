set terminal pngcairo font ",14"

set output "example1.png"

set xlabel "x"
set ylabel "u(x)"

fone(k,x) = (exp(k) - 1.0)/(exp(k) - exp(-k))*exp(-k*x) + (1.0 - exp(-k))/(exp(k) - exp(-k))*exp(k*x)

plot  fone(27.79,x) lw 2 title "analytical solution", "poisson_results.txt" every 5::0::100 u 1:2 w p pt 7 title "numerical (every 5th point)"


reset
set output "example2.png"

set xlabel "x"
set ylabel "u(x)"

ftwo(t,x) = cosh(t*(1.0-x))/cosh(t)

plot  ftwo(1.0,x) lw 2 title "analytical solution", "reaction_results.txt" every 5::0::100 u 1:2 w p pt 7 title "numerical (every 5th point)"
