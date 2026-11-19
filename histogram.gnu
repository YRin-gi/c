# 出力フォーマットの設定 (PNG画像として保存)
set terminal pngcairo size 400,300 font "sans,12"
set output "histogram_with_errors.png"

# グラフの見た目・タイトルの設定
set title "Histogram with Standard Errors"
set xlabel "Interval Center (x)"
set ylabel "Estimated Value (y)"
set grid y
set xrange [0:1]
set yrange [0:1.5]                      # Y軸のグリッド線を表示

# ヒストグラム（箱）のスタイル設定
set style fill solid 0 border # 箱の中を塗りつぶし、枠線を表示
box_width = 1.0 / 16                 # 箱の幅を指定（区間幅に合わせて調整）

# プロット処理 (boxes と yerrorbars を重ね描き)
plot "histogram.dat" using 2:3:(box_width) with boxes title "Estimate" lc rgb "#4C72B0", \
     "histogram.dat" using 2:3:4 with yerrorbars title "Standard Error" lc rgb "#111111" lw 1.5 ps 0
