#! /bin/bash

./print_cube.sh Army 10 PSNR_final >> cube_fpsnr_army_s10
./print_cube.sh Army 20 PSNR_final >> cube_fpsnr_army_s20
./print_cube.sh Army 40 PSNR_final >> cube_fpsnr_army_s40

./print_cube.sh DogDance 10 PSNR_final >> cube_fpsnr_dogdance_s10
./print_cube.sh DogDance 20 PSNR_final >> cube_fpsnr_dogdance_s20
./print_cube.sh DogDance 40 PSNR_final >> cube_fpsnr_dogdance_s40

./print_cube.sh Evergreen 10 PSNR_final >> cube_fpsnr_evergreen_s10
./print_cube.sh Evergreen 20 PSNR_final >> cube_fpsnr_evergreen_s20
./print_cube.sh Evergreen 40 PSNR_final >> cube_fpsnr_evergreen_s40

#./print_cube.sh Mequon 10 PSNR_final >> cube_fpsnr_mequon_s10
#./print_cube.sh Mequon 20 PSNR_final >> cube_fpsnr_mequon_s20
#./print_cube.sh Mequon 40 PSNR_final >> cube_fpsnr_mequon_s40

./print_cube.sh Walking 10 PSNR_final >> cube_fpsnr_walking_s10
./print_cube.sh Walking 20 PSNR_final >> cube_fpsnr_walking_s20
./print_cube.sh Walking 40 PSNR_final >> cube_fpsnr_walking_s40

./print_cube.sh Army 10 PSNR_basic >> cube_bpsnr_army_s10
./print_cube.sh Army 20 PSNR_basic >> cube_bpsnr_army_s20
./print_cube.sh Army 40 PSNR_basic >> cube_bpsnr_army_s40

./print_cube.sh DogDance 10 PSNR_basic >> cube_bpsnr_dogdance_s10
./print_cube.sh DogDance 20 PSNR_basic >> cube_bpsnr_dogdance_s20
./print_cube.sh DogDance 40 PSNR_basic >> cube_bpsnr_dogdance_s40

./print_cube.sh Evergreen 10 PSNR_basic >> cube_bpsnr_evergreen_s10
./print_cube.sh Evergreen 20 PSNR_basic >> cube_bpsnr_evergreen_s20
./print_cube.sh Evergreen 40 PSNR_basic >> cube_bpsnr_evergreen_s40

#./print_cube.sh Mequon 10 PSNR_basic >> cube_bpsnr_mequon_s10
#./print_cube.sh Mequon 20 PSNR_basic >> cube_bpsnr_mequon_s20
#./print_cube.sh Mequon 40 PSNR_basic >> cube_bpsnr_mequon_s40

./print_cube.sh Walking 10 PSNR_basic >> cube_bpsnr_walking_s10
./print_cube.sh Walking 20 PSNR_basic >> cube_bpsnr_walking_s20
./print_cube.sh Walking 40 PSNR_basic >> cube_bpsnr_walking_s40

./print_cube.sh Army 10 time >> cube_time_army_s10
./print_cube.sh Army 20 time >> cube_time_army_s20
./print_cube.sh Army 40 time >> cube_time_army_s40

./print_cube.sh DogDance 10 time >> cube_time_dogdance_s10
./print_cube.sh DogDance 20 time >> cube_time_dogdance_s20
./print_cube.sh DogDance 40 time >> cube_time_dogdance_s40

./print_cube.sh Evergreen 10 time >> cube_time_evergreen_s10
./print_cube.sh Evergreen 20 time >> cube_time_evergreen_s20
./print_cube.sh Evergreen 40 time >> cube_time_evergreen_s40

#./print_cube.sh Mequon 10 time >> cube_time_mequon_s10
#./print_cube.sh Mequon 20 time >> cube_time_mequon_s20
#./print_cube.sh Mequon 40 time >> cube_time_mequon_s40

./print_cube.sh Walking 10 time >> cube_time_walking_s10
./print_cube.sh Walking 20 time >> cube_time_walking_s20
./print_cube.sh Walking 40 time >> cube_time_walking_s40
