#! /bin/bash
# Prints a table with the results of a sequence given as input
#
# Tables are as follows:
#
# rank \ n_sim
#    040 080 120 160 375
# 04
# 08
# 12
# 16
# 20
#
# Usage
# ./print_table.sh SEQUENCE SIGMA [MEASURE]
#
# SEQUENCE can be Army, DogDance, Evergreen, Mequon, Walking
# SIGMA    can be 10, 20, 40
# MEASURE  can be [PSNR_final], PSNR_basic, RMSE_final, RMSE_basic, time

measure=${3:-PSNR_final}

t_11=`cat $1_s$2_2r04_2np040/measures | grep $measure | sed "s/-$measure = //"` 
t_12=`cat $1_s$2_2r04_2np080/measures | grep $measure | sed "s/-$measure = //"` 
t_13=`cat $1_s$2_2r04_2np120/measures | grep $measure | sed "s/-$measure = //"` 
t_14=`cat $1_s$2_2r04_2np160/measures | grep $measure | sed "s/-$measure = //"` 
t_15=`cat $1_s$2_2r04_2np375/measures | grep $measure | sed "s/-$measure = //"` 
t_21=`cat $1_s$2_2r08_2np040/measures | grep $measure | sed "s/-$measure = //"` 
t_22=`cat $1_s$2_2r08_2np080/measures | grep $measure | sed "s/-$measure = //"` 
t_23=`cat $1_s$2_2r08_2np120/measures | grep $measure | sed "s/-$measure = //"` 
t_24=`cat $1_s$2_2r08_2np160/measures | grep $measure | sed "s/-$measure = //"` 
t_25=`cat $1_s$2_2r08_2np375/measures | grep $measure | sed "s/-$measure = //"` 
t_31=`cat $1_s$2_2r12_2np040/measures | grep $measure | sed "s/-$measure = //"` 
t_32=`cat $1_s$2_2r12_2np080/measures | grep $measure | sed "s/-$measure = //"` 
t_33=`cat $1_s$2_2r12_2np120/measures | grep $measure | sed "s/-$measure = //"` 
t_34=`cat $1_s$2_2r12_2np160/measures | grep $measure | sed "s/-$measure = //"` 
t_35=`cat $1_s$2_2r12_2np375/measures | grep $measure | sed "s/-$measure = //"` 
t_41=`cat $1_s$2_2r16_2np040/measures | grep $measure | sed "s/-$measure = //"` 
t_42=`cat $1_s$2_2r16_2np080/measures | grep $measure | sed "s/-$measure = //"` 
t_43=`cat $1_s$2_2r16_2np120/measures | grep $measure | sed "s/-$measure = //"` 
t_44=`cat $1_s$2_2r16_2np160/measures | grep $measure | sed "s/-$measure = //"` 
t_45=`cat $1_s$2_2r16_2np375/measures | grep $measure | sed "s/-$measure = //"` 
t_51=`cat $1_s$2_2r20_2np040/measures | grep $measure | sed "s/-$measure = //"` 
t_52=`cat $1_s$2_2r20_2np080/measures | grep $measure | sed "s/-$measure = //"` 
t_53=`cat $1_s$2_2r20_2np120/measures | grep $measure | sed "s/-$measure = //"` 
t_54=`cat $1_s$2_2r20_2np160/measures | grep $measure | sed "s/-$measure = //"` 
t_55=`cat $1_s$2_2r20_2np375/measures | grep $measure | sed "s/-$measure = //"` 

t_11=`echo "scale=2; (10^2*$t_11 + 0.5) / 10^2" | bc`
t_12=`echo "scale=2; (10^2*$t_12 + 0.5) / 10^2" | bc`
t_13=`echo "scale=2; (10^2*$t_13 + 0.5) / 10^2" | bc`
t_14=`echo "scale=2; (10^2*$t_14 + 0.5) / 10^2" | bc`
t_15=`echo "scale=2; (10^2*$t_15 + 0.5) / 10^2" | bc`
t_21=`echo "scale=2; (10^2*$t_21 + 0.5) / 10^2" | bc`
t_22=`echo "scale=2; (10^2*$t_22 + 0.5) / 10^2" | bc`
t_23=`echo "scale=2; (10^2*$t_23 + 0.5) / 10^2" | bc`
t_24=`echo "scale=2; (10^2*$t_24 + 0.5) / 10^2" | bc`
t_25=`echo "scale=2; (10^2*$t_25 + 0.5) / 10^2" | bc`
t_31=`echo "scale=2; (10^2*$t_31 + 0.5) / 10^2" | bc`
t_32=`echo "scale=2; (10^2*$t_32 + 0.5) / 10^2" | bc`
t_33=`echo "scale=2; (10^2*$t_33 + 0.5) / 10^2" | bc`
t_34=`echo "scale=2; (10^2*$t_34 + 0.5) / 10^2" | bc`
t_35=`echo "scale=2; (10^2*$t_35 + 0.5) / 10^2" | bc`
t_41=`echo "scale=2; (10^2*$t_41 + 0.5) / 10^2" | bc`
t_42=`echo "scale=2; (10^2*$t_42 + 0.5) / 10^2" | bc`
t_43=`echo "scale=2; (10^2*$t_43 + 0.5) / 10^2" | bc`
t_44=`echo "scale=2; (10^2*$t_44 + 0.5) / 10^2" | bc`
t_45=`echo "scale=2; (10^2*$t_45 + 0.5) / 10^2" | bc`
t_51=`echo "scale=2; (10^2*$t_51 + 0.5) / 10^2" | bc`
t_52=`echo "scale=2; (10^2*$t_52 + 0.5) / 10^2" | bc`
t_53=`echo "scale=2; (10^2*$t_53 + 0.5) / 10^2" | bc`
t_54=`echo "scale=2; (10^2*$t_54 + 0.5) / 10^2" | bc`
t_55=`echo "scale=2; (10^2*$t_55 + 0.5) / 10^2" | bc`

printf "%s\t%s\t%s\t%s\t%s\n" $t_11 $t_12 $t_13 $t_14 $t_15
printf "%s\t%s\t%s\t%s\t%s\n" $t_21 $t_22 $t_23 $t_24 $t_25
printf "%s\t%s\t%s\t%s\t%s\n" $t_31 $t_32 $t_33 $t_34 $t_35
printf "%s\t%s\t%s\t%s\t%s\n" $t_41 $t_42 $t_43 $t_44 $t_45
printf "%s\t%s\t%s\t%s\t%s\n" $t_51 $t_52 $t_53 $t_54 $t_55

