# 出力フォーマットの設定 (PNG画像として保存)
set terminal pngcairo size 2000,1500 font "sans,12"
set output "2-2.png"

# グラフの見た目・タイトルの設定
set title "~e^-x/2"
set xlabel "x"
set ylabel "y"
set grid y
set autoscale         # Y軸のグリッド線を表示

# ヒストグラム（箱）のスタイル設定
set style fill solid 0 border # 箱の中を塗りつぶし、枠線を表示
box_width = 10.0 / 20                 # 箱の幅を指定（区間幅に合わせて調整）

# プロット処理 (boxes と yerrorbars を重ね描き)
plot "2-2.dat" using 2:3:(box_width) with boxes title "Estimate" lc rgb "#4C72B0", \
     "2-2.dat" using 2:3:4 with yerrorbars title "Standard Error" lc rgb "#111111" lw 1.5 ps 0
