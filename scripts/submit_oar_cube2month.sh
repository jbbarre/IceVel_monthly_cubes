#!/bin/bash
cd /summer/ice_speed/surface_flow_velocity/ANTARCTICA/MOSAIC/MONTHLY_MAPS/ASE_jbreproc

eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda activate cubes_py37
export MKL_NUM_THREADS=2

cpu=32
while read p; do
    i=$(echo "$p" | grep -oP 'x\K[0-9]+')
    j=$(echo "$p" | grep -oP 'y\K[0-9]+')
    echo $i $j
    python $HOME/ST_RELEASE/MOSAIC/PYTHON/submit_oar_cube2month.py $i $j > output$i-$j &
    cpu=`ps -fe | grep mougin | grep submit_oar_cube2month.py | wc -l`
    echo "# cpus used : " $cpu
    while [ $cpu -ge 16 ]
    do
        cpu=`ps -fe | grep mougin | grep submit_oar_cube2month.py | wc -l`
        echo "# cpus used : " $cpu
        if [ $cpu -ge 16 ]
        then
            sleep 3m
        fi
    done
done < /summer/ice_speed/surface_flow_velocity/ANTARCTICA/MOSAIC/MONTHLY_MAPS/ASE_jbreproc/$1

#wait until finished
cpu=`ps -fe | grep mougin | grep submit_oar_cube2month.py | wc -l`
while [ $cpu -ge 1 ]
    do
        cpu=`ps -fe | grep mougin | grep submit_oar_cube2month.py | wc -l`
        echo "# cpus used : " $cpu
        if [ $cpu -ge 1 ]
        then
            sleep 3m
        fi
done


