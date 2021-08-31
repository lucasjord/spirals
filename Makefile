#
# v1 : SDW First Makefile for the fortran code

FF=gfortran

all:
        ${FF} -o fit_geoblocks_tropos fit_geoblocks_tropos.f
        ${FF} -o fit_geoblocks_w_baselines fit_geoblocks_w_baselines.f
        ${FF} -o fit_parallax_multi fit_parallax_multi.f
        ${FF} -o  multiview  multiview.f

clean:
        rm fit_geoblocks_tropos fit_geoblocks_w_baselines fit_parallax_multi multiview_4x

install:
        cp fit_geoblocks_tropos ~/bin/fit_geoblocks_tropos
        cp fit_geoblocks_w_baselines ~/bin/fit_geoblocks_w_baselines
        cp fit_parallax_multi ~/bin/fit_parallax_multi
        cp multiview ~/bin/multiview

#End of the makefile
