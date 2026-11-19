reset

f(x)=a*x+b
a=1
b=1
fit f(x) 'pulse.txt' using 1:2 via a,b
print "V=",a,"E+",b
print "E=",1/a,"V+",b/a

# 出力形式を PNG に設定
set terminal pngcairo size 800,600
set output 'output.png'

set grid
set xlabel "E/keV"
set ylabel "V/V"
set title ""
plot f(x), 'pulse.txt' using 1:2
