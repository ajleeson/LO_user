#!/bin/bash

python ../../../LO/extract/box/extract_box.py -gtx cas7_t0_x4b -ro 5 -0 2014.01.01 -1 2014.12.31 -lt daily -job pugetsoundDO
python ../../../LO/extract/box/extract_box.py -gtx cas7_t0_x4b -ro 5 -0 2015.01.01 -1 2015.12.31 -lt daily -job pugetsoundDO
python ../../../LO/extract/box/extract_box.py -gtx cas7_t0_x4b -ro 5 -0 2016.01.01 -1 2016.12.31 -lt daily -job pugetsoundDO
python ../../../LO/extract/box/extract_box.py -gtx cas7_t0_x4b -ro 5 -0 2017.01.01 -1 2017.12.31 -lt daily -job pugetsoundDO
python ../../../LO/extract/box/extract_box.py -gtx cas7_t0_x4b -ro 5 -0 2018.01.01 -1 2018.12.31 -lt daily -job pugetsoundDO
python ../../../LO/extract/box/extract_box.py -gtx cas7_t0_x4b -ro 5 -0 2019.01.01 -1 2019.12.31 -lt daily -job pugetsoundDO