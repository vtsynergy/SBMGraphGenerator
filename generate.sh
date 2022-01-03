#! /bin/bash
directory="--directory=/groups/synergy_lab/wanyef/synthetic_graphs/"

vertices=100000

echo "==== Changes in sparsity ===="
let totalgraphs=10
let i=1
for density in 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 ; do
  echo "==== producing graph ${i} / ${totalgraphs} density = ${density} ===="
  python single_graph_generator.py -n ${vertices} -c -1 -m 0.35 -d ${density} \
      -o 2 -s 5 ${directory}
  let i+=1
done

echo "==== Changes in distribution ===="
let totalgraphs=8
let i=1
for exponent in -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5
do
  echo "==== producing graph ${i} / ${totalgraphs} exponent = ${exponent} ===="
  python single_graph_generator.py -n ${vertices} -c -1 -m 0.35 -o 2 -s 5 -e ${exponent} ${directory}
  let i+=1
done

echo "==== Changes in community size ===="
let totalgraphs=12
let i=1
let size=32768  # power of 2
for comm in 4096 2048 1024 512 256 128 64 32 16 8 4 2  # 1000 2000 4000 8000 16000 32000 64000 128000 256000 512000 1024000
do
  echo "====producing graph ${i} / ${totalgraphs} communities = ${comm} ===="
  python single_graph_generator.py -n ${size} -c ${comm} -o 2 -s 5 ${directory}
  let i+=1
done

echo "==== Changes in community overlap ===="
let totalgraphs=11
let i=1
for overlap in 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0  # -3.5 -2.8 -2.1 -1.4 -0.7 0 0.7 1.4 2.1 2.8 3.5 
do
  echo "==== producing graph ${i} / ${totalgraphs} overlap = ${overlap} ===="
  python single_graph_generator.py -n ${vertices} -c -1 -m 0.35 -o ${overlap} -s 5 ${directory}
  let i+=1
done

echo "==== Changes in community size heterogeneity ===="
let totalgraphs=11
let i=1
for blockvar in 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0  # -3.5 -2.8 -2.1 -1.4 -0.7 0 0.7 1.4 2.1 2.8 3.5 
do
  echo "==== producing graph ${i} / ${totalgraphs} blockvar = ${blockvar} ===="
  python single_graph_generator.py -n ${vertices} -c -1 -m 0.35 -o 2 -s ${blockvar} ${directory}
  let i+=1
done

echo "==== Changes in graph size ===="
let totalgraphs=12
let i=1
for size in 1000 2000 4000 8000 16000 32000 64000 128000 256000 512000 1024000 2048000
do
  echo "====producing graph ${i} / ${totalgraphs} size = ${size} ===="
  python single_graph_generator.py -n ${size} -c -1 -m 0.35 -o 2 -s 5 ${directory}
  let i+=1
done

echo "==== Changes in supernode strength ===="
let mindegree=1
let totalgraphs=7
let i=1
for maxdegree in 0.0005 0.001 0.005 0.01 0.05 0.1 0.5
do
  echo "==== producing graph ${i} / ${totalgraphs} max = ${maxdegree} ===="
  python single_graph_generator.py -n ${vertices} -c -1 -m 0.35 -a ${maxdegree} -o 2 -s 5 ${directory}
  let i+=1
done

