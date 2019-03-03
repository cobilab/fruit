# !/usr/bin/env bash

### To activate, set to 1
DATASET_SYNTH=1
DATASET_ANCIENT=1
DATASET_HS=1
RUN_MUT_SYNTH=1
RUN_MUT_ANCIENT=1
RUN_HS=1

########## DATASET - SYNTHETIC ##########
if [[ $DATASET_SYNTH -eq 1 ]]; then
  if [[ ! -d ../dataset ]]; then mkdir -p ../dataset; fi
  if [[ ! -d ../dataset/synth ]]; then mkdir -p ../dataset/synth; fi
  cp ../prog/mutate ../dataset/synth
  cd ../dataset/synth

  < /dev/urandom tr -d -c 'A-T' | head -c 5000000 > syn

  for i in {1..50}; do
    ./mutate -s $i -r $i -i syn -o mut_$i
  done
fi

########## DATASET - ANCIENT ##########
if [[ $DATASET_ANCIENT -eq 1 ]]; then
  cd ..
  if [[ ! -d dataset ]]; then mkdir -p dataset; fi
  if [[ ! -d dataset/ancient ]]; then mkdir -p dataset/ancient; fi
  cd dataset/ancient
  cp ../../prog/goose-fasta2seq ../../prog/mutate  .

  wget --trust-server-names http://cdna.eva.mpg.de/neandertal/exomes/proteins/primary_target.transcripts.CDS_3.Altai.prot.fa
  mv primary_target.transcripts.CDS_3.Altai.prot.fa  A.fa
  
  wget --trust-server-names http://cdna.eva.mpg.de/neandertal/exomes/proteins/primary_target.transcripts.CDS_3.Sidron.prot.fa
  mv primary_target.transcripts.CDS_3.Sidron.prot.fa  S.fa

  wget --trust-server-names http://cdna.eva.mpg.de/neandertal/exomes/proteins/primary_target.transcripts.CDS_3.Vindija.prot.fa
  mv primary_target.transcripts.CDS_3.Vindija.prot.fa  V.fa

  ### Convert to Seq
  ./goose-fasta2seq < A.fa > A.seq
  ./goose-fasta2seq < S.fa > S.seq
  ./goose-fasta2seq < V.fa > V.seq

  ### Mutate
  for i in {0..50}; do
    ./mutate -s $i -r $i -i A.seq -o A_mut_$i
    ./mutate -s $i -r $i -i S.seq -o S_mut_$i
    ./mutate -s $i -r $i -i V.seq -o V_mut_$i
  done

  rm -f mutate goose-fasta2seq
  cd ../../script
fi

########## DATASET - HS ##########
if [[ $DATASET_HS -eq 1 ]]; then
  if [[ ! -d dataset ]]; then mkdir -p dataset; fi
  if [[ ! -d dataset/hs ]]; then mkdir -p dataset/hs; fi
  cd ../dataset/ancient
  cat A.seq S.seq V.seq > ../hs/ASV.seq
  cd ../..
  cp prog/goose-splitreads prog/goose-fasta2seq HS.mfasta  dataset/hs
  cd dataset/hs

  # ##todo: doesn't work
  # # wget --trust-server-names https://www.uniprot.org/uniprot/?query=*&format=fasta&force=true&fil=proteome:UP000005640%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22&compress=no
  # # mv uniprot-filtered-proteome%3AUP000005640+AND+reviewed%3Ayes+AND+organism%3A_Hom--.fasta  HS.fasta

  ./goose-splitreads < HS.mfasta

  for file in *.fasta; do
    ./goose-fasta2seq < $file > ${file%.*}.seq
  done

  ## Extract names
  rm -f names
  for file in *.seq; do
    echo -e "${file/.seq/}\t$(head -n1 $file)" >> names
  done
  sort -k1 -V names -o names
  cp names ../../result/hs

  rm -f goose-splitreads goose-fasta2seq HS.mfasta *.fasta
fi

########## RUN MUTATED SYNTHETIC ##########
if [[ $RUN_MUT_SYNTH -eq 1 ]]; then
  cd ..
  sudo chmod -R 777 .
  cp fruit-map fruit-filter fruit-visual dataset/synth
  cd dataset/synth

  MIN_MUTE=0
  MAX_MUTE=50
  MIN_K=1
  MAX_K=12

  ### Map
  target=""
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r syn -t $target -k $i -nth 8

    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      mv syn-mut_${j}.sng  ${i}-syn-mut_${j}.sng;
    done
  done

  ### Filter
  target=""
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      target="$target:${i}-syn-mut_${j}.sng"
    done
  done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-filter -i $target -w 1 -nth 8
  done

  echo -e "k\tmut\tname\tunqRat\tunqCnt\tnotUnqCnt\ttotCnt" > uniq
  echo -e "k\tunqRat\tunqCnt\ttotCnt" > stat
  for (( KMER=MIN_K; KMER<=MAX_K; KMER++ )); do
    echo "Calculating uniqueness ratios for k-mer=$KMER ..."
    rm -f uniq$KMER
    totUniqCount=0
    totSize=0

    for file in $KMER-*.pos; do
      {
      read watermark name size
      uniqCount=0
      while read begin end; do
        if [[ $begin != "" ]]; then
          uniqCount=$((uniqCount+end-begin+1))
        fi
      done
      uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
      notUniqCount=$((size-uniqCount))
      echo -e "$name\t$uniqRatio\t$uniqCount\t$notUniqCount\t$size" >> uniq$KMER
      totUniqCount=$((totUniqCount+uniqCount))
      totSize=$((totSize+size))
      } < $file
    done

    totUniqRatio=$(bc <<< "scale=5; $totUniqCount/$totSize")
    echo -e "$KMER\t$totUniqRatio\t$totUniqCount\t$totSize" >> stat

    ## Sort
    if [[ -f uniq$KMER ]]; then
      sort -k1,1V -o uniq$KMER uniq$KMER

      mut=0
      while IFS="\t" read all; do
        echo -e "$KMER\t$mut\t$all" >> uniq
        mut=$((mut+1))
      done < uniq$KMER
    fi

    rm -f uniq$KMER *.sng
    echo "Finished.";
  done

  if [[ ! -d ../../result/synth ]]; then mkdir -p ../../result/synth; fi
  cp uniq stat  ../../result/synth

  cd ../../script
fi

########## RUN MUTATED ANCIENT ##########
if [[ $RUN_MUT_ANCIENT -eq 1 ]]; then
  cd ..
  sudo chmod -R 777 .
  cp fruit-map fruit-filter fruit-visual  dataset/ancient
  cd dataset/ancient

  MIN_MUTE=0
  MAX_MUTE=50
  MIN_K=1
  MAX_K=12

  ### Altai
  ## Map
  target=""  
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:A_mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r A.seq -t $target -k $i -nth 8

    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      mv A.seq-A_mut_${j}.sng  ${i}-A.seq-A_mut_${j}.sng;
    done
  done
  
  ## Filter
  target=""
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      target="$target:${i}-A.seq-A_mut_${j}.sng"
    done
  done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-filter -i $target -w 1 -nth 8
  done

  echo -e "k\tmut\tname\tunqRat\tunqCnt\tnotUnqCnt\ttotCnt" > uniq-A
  echo -e "k\tunqRat\tunqCnt\ttotCnt" > stat-A
  for (( KMER=MIN_K; KMER<=MAX_K; KMER++ )); do
    echo "Calculating uniqueness ratios for k-mer=$KMER ..."
    rm -f uniq$KMER
    totUniqCount=0
    totSize=0

    for file in $KMER-A.seq-*.pos; do
      {
      read watermark name size
      uniqCount=0
      while read begin end; do
        if [[ $begin != "" ]]; then
          uniqCount=$((uniqCount+end-begin+1))
        fi
      done
      uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
      notUniqCount=$((size-uniqCount))
      echo -e "$name\t$uniqRatio\t$uniqCount\t$notUniqCount\t$size" >> uniq$KMER
      totUniqCount=$((totUniqCount+uniqCount))
      totSize=$((totSize+size))
      } < $file
    done

    totUniqRatio=$(bc <<< "scale=5; $totUniqCount/$totSize")
    echo -e "$KMER\t$totUniqRatio\t$totUniqCount\t$totSize" >> stat-A

    ## Sort
    if [[ -f uniq$KMER ]]; then
      sort -k1,1V -o uniq$KMER uniq$KMER

      mut=0
      while IFS="\t" read all; do
        echo -e "$KMER\t$mut\t$all" >> uniq-A
        mut=$((mut+1))
      done < uniq$KMER
    fi

    rm -f uniq$KMER *.sng
    echo "Finished.";
  done

  ### Sidron
  ## Map
  target=""
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:S_mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r S.seq -t $target -k $i -nth 8

    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      mv S.seq-S_mut_${j}.sng  ${i}-S.seq-S_mut_${j}.sng;
    done
  done
  
  ## Filter
  target=""
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      target="$target:${i}-S.seq-S_mut_${j}.sng"
    done
  done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-filter -i $target -w 1 -nth 8
  done

  echo -e "k\tmut\tname\tunqRat\tunqCnt\tnotUnqCnt\ttotCnt" > uniq-S
  echo -e "k\tunqRat\tunqCnt\ttotCnt" > stat-S
  for (( KMER=MIN_K; KMER<=MAX_K; KMER++ )); do
    echo "Calculating uniqueness ratios for k-mer=$KMER ..."
    rm -f uniq$KMER
    totUniqCount=0
    totSize=0

    for file in $KMER-S.seq-*.pos; do
      {
      read watermark name size
      uniqCount=0
      while read begin end; do
        if [[ $begin != "" ]]; then
          uniqCount=$((uniqCount+end-begin+1))
        fi
      done
      uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
      notUniqCount=$((size-uniqCount))
      echo -e "$name\t$uniqRatio\t$uniqCount\t$notUniqCount\t$size" >> uniq$KMER
      totUniqCount=$((totUniqCount+uniqCount))
      totSize=$((totSize+size))
      } < $file
    done

    totUniqRatio=$(bc <<< "scale=5; $totUniqCount/$totSize")
    echo -e "$KMER\t$totUniqRatio\t$totUniqCount\t$totSize" >> stat-S

    ## Sort
    if [[ -f uniq$KMER ]]; then
      sort -k1,1V -o uniq$KMER uniq$KMER

      mut=0
      while IFS="\t" read all; do
        echo -e "$KMER\t$mut\t$all" >> uniq-S
        mut=$((mut+1))
      done < uniq$KMER
    fi

    rm -f uniq$KMER *.sng
    echo "Finished."
  done

  ### Vindija
  ## Map
  target=""
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:V_mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r V.seq -t $target -k $i -nth 8

    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      mv V.seq-V_mut_${j}.sng  ${i}-V.seq-V_mut_${j}.sng;
    done
  done
  
  ## Filter
  target=""
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      target="$target:${i}-V.seq-V_mut_${j}.sng"
    done
  done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-filter -i $target -w 1 -nth 8
  done

  echo -e "k\tmut\tname\tunqRat\tunqCnt\tnotUnqCnt\ttotCnt" > uniq-V
  echo -e "k\tunqRat\tunqCnt\ttotCnt" > stat-V
  for (( KMER=MIN_K; KMER<=MAX_K; KMER++ )); do
    echo "Calculating uniqueness ratios for k-mer=$KMER ..."
    rm -f uniq$KMER
    totUniqCount=0
    totSize=0

    for file in $KMER-V.seq-*.pos; do
      {
      read watermark name size
      uniqCount=0
      while read begin end; do
        if [[ $begin != "" ]]; then
          uniqCount=$((uniqCount+end-begin+1))
        fi
      done
      uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
      notUniqCount=$((size-uniqCount))
      echo -e "$name\t$uniqRatio\t$uniqCount\t$notUniqCount\t$size" >> uniq$KMER
      totUniqCount=$((totUniqCount+uniqCount))
      totSize=$((totSize+size))
      } < $file
    done

    totUniqRatio=$(bc <<< "scale=5; $totUniqCount/$totSize")
    echo -e "$KMER\t$totUniqRatio\t$totUniqCount\t$totSize" >> stat-V

    ## Sort
    if [[ -f uniq$KMER ]]; then
      sort -k1,1V -o uniq$KMER uniq$KMER

      mut=0
      while IFS="\t" read all; do
        echo -e "$KMER\t$mut\t$all" >> uniq-V
        mut=$((mut+1))
      done < uniq$KMER
    fi

    rm -f uniq$KMER *.sng
    echo "Finished."
  done

  ### Merge
  echo -e "dataset\tk\tmut\tname\tunqRat\tunqCnt\tnotUnqCnt\ttotCnt" > uniq
  sed 1,1d uniq-A | while IFS="\t" read all; do 
    echo -e "Altai\t$all" >> uniq; 
  done
  sed 1,1d uniq-S | while read all; do 
    echo -e "Sidron\t$all" >> uniq; 
  done
  sed 1,1d uniq-V | while read all; do 
    echo -e "Vindija\t$all" >> uniq; 
  done

  echo -e "dataset\tk\tunqRat\tunqCnt\ttotCnt" > stat
  sed 1,1d stat-A | while IFS="\t" read all; do 
    echo -e "Altai\t$all" >> stat; 
  done
  sed 1,1d stat-S | while read all; do 
    echo -e "Sidron\t$all" >> stat; 
  done
  sed 1,1d stat-V | while read all; do 
    echo -e "Vindija\t$all" >> stat; 
  done

  rm -f uniq-A uniq-S uniq-V stat-A stat-S stat-V

  if [[ ! -d ../../result/ancient ]]; then mkdir -p ../../result/ancient; fi
  cp uniq stat  ../../result/ancient

  rm -f fruit-map fruit-filter fruit-visual
  cd ../../script
fi

########## RUN HS ##########
if [[ $RUN_HS -eq 1 ]]; then
  cd ..
  sudo chmod -R 777 .
  cp fruit-map fruit-filter fruit-visual  dataset/hs
  cd dataset/hs

  MIN=1
  MAX=20412
  MIN_K=1
  MAX_K=12
  STEP=8000
  TOP=20

  for (( KMER=$MIN_K; KMER<=$MAX_K; KMER++ )); do
    ### Map
    for (( i=MIN; i<=MAX; i+=STEP )); do
      target=""
      for (( j=i; j<=MAX && j<i+STEP; j++ )); do 
        target="$target:out$j.seq"; 
      done
      target=${target:1}

      ./fruit-map -r ASV.seq -t $target -k $KMER -bp 0.00001 -nth 8
    done

    for (( i=MIN; i<=MAX; i++ )); do
      mv ASV.seq-out${i}.seq.sng $KMER-out${i}.sng;
    done

    ### Filter
    for (( i=MIN; i<=MAX; i+=STEP )); do
      target=""
      for (( j=i; j<=MAX && j<i+STEP; j++ )); do 
        target="$target:$KMER-out$j.sng";
      done
      target=${target:1}

      ./fruit-filter -i $target -w 1 -nth 8
    done

    echo "Calculating uniqueness ratios ..."
    rm -f uniq$KMER
    totUniqCount=0
    totSize=0

    for file in $KMER-*.pos; do
      {
      read watermark name size
      uniqCount=0
      while read begin end; do
        if [[ $begin != "" ]]; then
          uniqCount=$((uniqCount+end-begin+1))
        fi
      done
      if [[ $uniqCount != 0 ]]; then
        uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
        notUniqCount=$((size-uniqCount))
        echo -e "$name\t$uniqRatio\t$uniqCount\t$notUniqCount\t$size" >> uniq$KMER
      fi
      totUniqCount=$((totUniqCount+uniqCount))
      totSize=$((totSize+size))
      } < $file
    done

    totUniqRatio=$(bc <<< "scale=5; $totUniqCount/$totSize")  
    echo -e "$KMER\t$totUniqRatio\t$totUniqCount\t$totSize" >> stat$KMER

    echo "Finished";
  done
  
  ## Sort & merge
  echo -e "k\tname\tunqRat\tunqCnt\tnotUnqCnt\ttotCnt" > uniq
  echo -e "k\tunqRat\tunqCnt\ttotCnt" > stat
  for (( KMER=$MIN_K; KMER<=$MAX_K; KMER++ )); do
    if [[ -f uniq$KMER ]]; then
      sort -k1,1V -o uniq$KMER uniq$KMER

      while IFS="\t" read all; do 
        echo -e "$KMER\t$all" >> uniq
      done < uniq$KMER

      rm -f uniq$KMER
    fi

    if [[ -f stat$KMER ]]; then
      cat stat$KMER >> stat
      rm -f stat$KMER
    fi
  done

  if [[ ! -d ../../result/hs ]]; then mkdir -p ../../result/hs; fi
  cp uniq stat  ../../result/hs

  ## Extract top and bottom for visualization
  head -n1 uniq > top
  awk '{ if ($1 == 7) { print } }' uniq >> vis
  sort -k3,3r -k2,2V vis | head -n $TOP >> top
  sort -k3,3 -k2,2V vis | head -n $TOP >> top
  rm -f vis

  ### Visual
  target=""
  target=$(awk -F"\t" '{printf ":"$2".pos"}' top);
  target=${target:1}
  # target=${target:10}
  ./fruit-visual -i $target -o top.svg

  ./fruit-visual -t 20 -am 280 -th 0 -b 70 -rh 4 -w 11 -s 20 -o top.svg -i 7-out4754.pos:7-out645.pos:7-out19770.pos:7-out16394.pos:7-out6067.pos:7-out20045.pos:7-out19013.pos:7-out19655.pos:7-out1070.pos:7-out4065.pos  -sl "Putative uncharacterized protein encoded by LINC00303":"Protein PAPPAS":"Ataxin-8":"Putative uncharacterized protein encoded by LINC01545":"60S ribosomal protein L26":"Uncharacterized protein C12orf71":"Interferon-induced transmembrane protein 3":"Uncharacterized protein C15orf65":"NADH-ubiquinone oxidoreductase chain 3":"Cytochrome c oxidase subunit 2"

  if [[ ! -d ../../result/hs ]]; then mkdir -p ../../result/hs; fi
  cp *.svg  ../../result/hs

  cd ../../script
fi