#!/bin/bash

pysph run taylor_green --kernel=CubicSpline --hdx=0.5 --nx=25 --tf=2.0 --openmp -d ./Cubicspline/case1
pysph run taylor_green --kernel=CubicSpline --hdx=0.5 --nx=50 --tf=2.0 --openmp -d ./Cubicspline/case2
pysph run taylor_green --kernel=CubicSpline --hdx=0.5 --nx=100 --tf=2.0 --openmp -d ./Cubicspline/case3
pysph run taylor_green --kernel=CubicSpline --hdx=1.0 --nx=25  --tf=2.0 --openmp -d ./Cubicspline/case4
pysph run taylor_green --kernel=CubicSpline --hdx=1.0 --nx=50 --tf=2.0 --openmp -d  ./Cubicspline/case5
pysph run taylor_green --kernel=CubicSpline --hdx=1.0 --nx=100 --tf=2.0 --openmp -d ./Cubicspline/case6
pysph run taylor_green --kernel=CubicSpline --hdx=2.0 --nx=25 --tf=2.0 --openmp -d  ./Cubicspline/case7
pysph run taylor_green --kernel=CubicSpline --hdx=2.0 --nx=50 --tf=2.0 --openmp -d  ./Cubicspline/case8
pysph run taylor_green --kernel=CubicSpline --hdx=2.0 --nx=100 --tf=2.0 --openmp -d ./Cubicspline/case9
pysph run taylor_green --kernel=Gaussian --hdx=0.5 --nx=25 --tf=2.0 --openmp -d ./gaussian/case1
pysph run taylor_green --kernel=Gaussian --hdx=0.5 --nx=50 --tf=2.0 --openmp -d ./gaussian/case2
pysph run taylor_green --kernel=Gaussian --hdx=0.5 --nx=100 --tf=2.0 --openmp -d ./gaussian/case3
pysph run taylor_green --kernel=Gaussian --hdx=1.0 --nx=25 --tf=2.0 --openmp -d ./gaussian/case4
pysph run taylor_green --kernel=Gaussian --hdx=1.0 --nx=50 --tf=2.0 --openmp -d ./gaussian/case5
pysph run taylor_green --kernel=Gaussian --hdx=1.0 --nx=100 --tf=2.0 --openmp -d ./gaussian/case6
pysph run taylor_green --kernel=Gaussian --hdx=2.0 --nx=25 --tf=2.0 --openmp -d ./gaussian/case7
pysph run taylor_green --kernel=Gaussian --hdx=2.0 --nx=50 --tf=2.0 --openmp -d ./gaussian/case8
pysph run taylor_green --kernel=Gaussian --hdx=2.0 --nx=100 --tf=2.0 --openmp -d ./gaussian/case9
pysph run taylor_green --kernel=QuinticSpline --hdx=0.5 --nx=25 --tf=2.0 --openmp -d ./quintic/case1
pysph run taylor_green --kernel=QuinticSpline --hdx=0.5 --nx=50 --tf=2.0 --openmp -d ./quintic/case2
pysph run taylor_green --kernel=QuinticSpline --hdx=0.5 --nx=100 --tf=2.0 --openmp -d ./quintic/case3
pysph run taylor_green --kernel=QuinticSpline --hdx=1.0 --nx=25 --tf=2.0 --openmp -d ./quintic/case4
pysph run taylor_green --kernel=QuinticSpline --hdx=1.0 --nx=50 --tf=2.0 --openmp -d ./quintic/case5
pysph run taylor_green --kernel=QuinticSpline --hdx=1.0 --nx=100 --tf=2.0 --openmp -d ./quintic/case6
pysph run taylor_green --kernel=QuinticSpline --hdx=2.0 --nx=25 --tf=2.0 --openmp -d ./quintic/case7
pysph run taylor_green --kernel=QuinticSpline --hdx=2.0 --nx=50 --tf=2.0 --openmp -d ./quintic/case8
pysph run taylor_green --kernel=QuinticSpline --hdx=2.0 --nx=100 --tf=2.0 --openmp -d ./quintic/case9
pysph run taylor_green --kernel=WendlandQuintic --hdx=0.5 --nx=25 --tf=2.0 --openmp -d ./WendlandQuintic/case1
pysph run taylor_green --kernel=WendlandQuintic --hdx=0.5 --nx=50 --tf=2.0 --openmp -d ./WendlandQuintic/case2
pysph run taylor_green --kernel=WendlandQuintic --hdx=0.5 --nx=100 --tf=2.0 --openmp -d ./WendlandQuintic/case3
pysph run taylor_green --kernel=WendlandQuintic --hdx=1.0 --nx=25  --tf=2.0 --openmp -d ./WendlandQuintic/case4
pysph run taylor_green --kernel=WendlandQuintic --hdx=1.0 --nx=50 --tf=2.0 --openmp -d ./WendlandQuintic/case5
pysph run taylor_green --kernel=WendlandQuintic --hdx=1.0 --nx=100 --tf=2.0 --openmp -d ./WendlandQuintic/case6
pysph run taylor_green --kernel=WendlandQuintic --hdx=2.0 --nx=25 --tf=2.0 --openmp -d ./WendlandQuintic/case7
pysph run taylor_green --kernel=WendlandQuintic --hdx=2.0 --nx=50 --tf=2.0 --openmp -d ./WendlandQuintic/case8
pysph run taylor_green --kernel=WendlandQuintic --hdx=2.0 --nx=100 --tf=2.0 --openmp -d ./WendlandQuintic/case9

python a10-130010058.py
pdflatex a10-130010058.tex
pdflatex a10-130010058.tex
rm a10-130010058.aux
rm a10-130010058.log
rm a10-130010058.lol
rm a10-130010058.out
rm a10-130010058.toc
rm -rf ./CubicSpline
rm -rf ./gaussian
rm -rf ./quintic
rm -rf ./WendlandQuintic
rm *.png
