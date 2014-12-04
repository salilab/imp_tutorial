#!/bin/bash
# Simple script to plot data from IMP stat files using gnuplot
#
# usage:
#
# ./plot_stat.sh -sg -i STATFILE -y YCOLUMN [-x XCOLUMN] [-m POINTS] [-plot] [-o OUTPUTFILE] [-b BEGIN]
# 
# -i | input stat file name
# -x | column number with Y data values OR column header string [no x input means only plot y]
# -y | column number with Y data values OR column header string
# -m | method of plotting.  POINTS, LINES or LINESPOINTS
# -s | suppress showing plot
# -g | save gnuplot file
# -b | BEGIN from this frame (in case you want to skip an equilibration phase)
# -o | saves plot to OUTPUTFILE. default is no output
# -h | prints help text to screen
#
#


function get_col_field(){

clm=$1

if [[ $clm =~ ^-?[0-9]+$ ]]; then
    f=`echo $header | awk -v x=$clm -F"," '{print $x}' | awk '{print $NF}'`
else
    c=`echo $header | awk -v x=$clm -F"," '{
        for (i=1; i<=NF; i++)
            if (match($i, x) != 0) {
                print i
                exit
            }
        }'`
    clm=$c
    f=`echo $header | awk -v x=$clm -F"," '{print $x}' | awk '{print $NF}'`
fi

echo $clm $f


}


helptext=" # \n
# plot_stat.sh plots two columns from an IMP produced stat file. \n
# \n
# usage: \n
# \n
# ./plot_stat.sh -i STATFILE -x XCOLUMN -y YCOLUMN [-m POINTS] [-plot] [-o OUTPUTFILE] \n
# \n
# -i | input stat file name \n
# -x | column number with Y data values OR column header string \n
# -y | column number with Y data values OR column header string \n
# -m | method of plotting.  POINTS, LINES or LINESPOINTS \n
# -s | suppress showing plot \n
# -o | saves plot to OUTPUTFILE default is no output \n
# -h | prints this help text to screen \n
#\n
#\n"

# defaults
method="points"
output=false
x_input=false
show_plot=true
x_col=0
y_col=0
begin=1


# Parse options
while getopts ":hsgi:x:y:m:b:" flag; do
    case $flag in
       h)  echo -e $helptext; exit;;
       i)  stat_file=$OPTARG;;
       x)  x_input=true; x_col=$OPTARG;;
       y)  y_col=$OPTARG;;
       m)  method=$OPTARG;;
       b)  begin=$OPTARG;;
       s)  show_plot=false;;
       g)  save_gnuplot_file=true;;
       o)  save_plot=true; outputfile=$OPTARG;;
       *)  echo "Invalid Option ";;
    esac
done

header=`head -n1 $stat_file | cut -d "{" -f2 | cut -d "}" -f1`

#--------------------------
# Get data from stat file
#--------------------------

# If x_col or y_col are integers, grab that column header
# If they are strings, find the column ID
n_lines=`wc -l $stat_file | awk '{print $1}'`
n_lines=$(($n_lines-1))

if $x_input; then
    results=`get_col_field $x_col`
    x_col=`echo $results | awk '{print $1}'`
    x_field=`echo $results | awk '{print $2}' | tr "'" ' ' | sed 's/ //g'`
    echo "x values from column" $x_col":  " $x_field
    x_vals=`awk -v x=$x_col -F"," '{print $x}' $stat_file | awk '{print $NF}' | sed "s/'//g" | tail -n+2`
else  
    echo "x values is frame number"
    x_field="Frames"
    x_vals=`seq 1 $n_lines`
fi

results=`get_col_field $y_col`
y_col=`echo $results | awk '{print $1}'`
y_field=`echo $results | awk '{print $2}' | tr "'" ' '| sed 's/ //g'`
y_vals=`awk -v y=$y_col -F"," '{print $y}' $stat_file | awk '{print $NF}' | sed "s/'//g" | tail -n+2`
echo "y values from column" $y_col":  " $y_field

paste <(echo "$x_vals") <(echo "$y_vals") --delimiters " " > temp_stat.dat

# Create gnuplot file

gnuplot_file="test.gnuplot"
outfile=$x_field-$y_field.eps
if $show_plot; then
 echo "set terminal dumb" > $gnuplot_file
else
 echo "set terminal postscript enhanced" > $gnuplot_file
 echo "set output \"$outfile\"" >> $gnuplot_file
 echo "plot saved to " $outfile
fi
echo "unset key" >> $gnuplot_file
echo "set xlabel \"$x_field\"" >> $gnuplot_file
echo "set ylabel \"$y_field\"" >> $gnuplot_file
echo "plot \"temp_stat.dat\" every ::$begin::$n_lines u 1:2 w $method">> $gnuplot_file

gnuplot --persist $gnuplot_file

#rm temp_stat.dat
#rm $gnuplot_file




