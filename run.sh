# !/usr/bin/env bash

DATASET_SYNTH=1
DATASET_ANCIENT=1
DATASET_HS=1
RUN_MUT_SYNTH=1
RUN_MUT_ANCIENT=1
RUN_HS=1

########## DATASET - SYNTHETIC - 250 MB ##########
if [[ $DATASET_SYNTH -eq 1 ]]; then
  cp mutate dataset/synth
  cd dataset/synth

  # < /dev/urandom tr -d -c 'A-T' | head -c 5000000 > syn

  echo "Mutating synthetic dataset ..."
  for i in {1..50}; do
    ./mutate -s $i -r $i -i syn -o mut_$i
  done
  echo "finished."

  rm -f mutate
  cd ../..
fi

########## DATASET - ANCIENT - 3.6 GB ##########
if [[ $DATASET_ANCIENT -eq 1 ]]; then
  cp mutate dataset/ancient
  cd dataset/ancient

  echo "Mutating ancient human dataset ..."
  for i in {1..50}; do
    ./mutate -s $i -r $i -i A.fa -o A_mut_$i
    ./mutate -s $i -r $i -i S.fa -o S_mut_$i
    ./mutate -s $i -r $i -i V.fa -o V_mut_$i
  done
  echo "finished."

  rm -f mutate
  cd ../..
fi

########## DATASET - HS - 26 MB ##########
if [[ $DATASET_HS -eq 1 ]]; then
  cp goose-splitreads dataset/hs
  cd dataset/hs

  echo "Splitting modern human dataset ..."
  ./goose-splitreads < HS.mfasta
  echo "finished."

  echo "Renaming ..."
  for file in *.fasta; do
    mv "$file" "${file/.fasta/.fa}"
  done
  echo "finished."

  echo "Extracting names ..."
  rm -f names
  for file in *.fa; do
    echo -e "${file/.fa/}\t$(head -n1 $file)" >> names
  done
  sort -k1 -V names -o names
  echo "finished."

  if [[ ! -d ../../result ]]; then mkdir -p ../../result; fi
  cp names ../../result

  rm -f goose-splitreads
  cd ../..
fi

########## RUN MUTATED SYNTHETIC ##########
if [[ $RUN_MUT_SYNTH -eq 1 ]]; then
  cp fruit-map fruit-filter fruit-visual  dataset/synth
  cd dataset/synth

  MIN_MUTE=1
  MAX_MUTE=50
  MIN_K=6
  MAX_K=12

  echo "Running FRUIT on mutated synthetic dataset syn ..."
  ### Map
  target=""
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r syn -t $target -k $i

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
    ./fruit-filter -i $target -w 5
  done

  ## Remove empty files
  for file in *.pos; do
    numLines=$(wc -l < $file)
    if [[ $numLines -eq 1 ]]; then rm -f $file; fi
  done

  echo "Calculating uniqueness ratio ..."
  uniqRatio=0
  for file in *.pos; do
    {
    read watermark name size
    uniqCount=0
    while read begin end; do
      uniqCount=$((uniqCount+end-begin+1))
    done
    uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
    printf '%.5f' "$uniqRatio" > $file.unr
    } < $file
  done

  ## Save result
  echo -e "#mut\tuniqRat\tk" > "result-synth"
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      printf '%d\t%.5f\t%d\n' "$j" "$(< ${i}-syn-mut_${j}.pos.unr)" "$i" >> "result-synth"
    done
  done
  echo "Finished."

  rm -f *.unr *.sng

  if [[ ! -d ../../result ]]; then mkdir -p ../../result; fi
  cp result-synth  ../../result;

  rm -f fruit-map fruit-filter fruit-visual
  cd ../..
fi

########## RUN MUTATED ANCIENT ##########
if [[ $RUN_MUT_ANCIENT -eq 1 ]]; then
  cp fruit-map fruit-filter fruit-visual dataset/ancient
  cd dataset/ancient

  MIN_MUTE=1
  MAX_MUTE=50
  MIN_K=10
  MAX_K=12
  
  ### Altai
  echo "Running FRUIT on mutated ancient dataset Altai ..."

  target=""  
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:A_mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r A.fa -t $target -k $i

    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      mv A.fa-A_mut_${j}.sng  ${i}-A.fa-A_mut_${j}.sng;
    done
  done
  
  ## Filter
  target=""
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      target="$target:${i}-A.fa-A_mut_${j}.sng"
    done
  done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-filter -i $target -w 5
  done

  # Remove empty files
  for file in *A.fa*.pos; do
    numLines=$(wc -l < $file)
    if [[ $numLines -eq 1 ]]; then rm -f $file; fi
  done

  echo "Calculating uniqueness ratio ..."
  uniqRatio=0
  for file in *A.fa*.pos; do
    {
    read watermark name size
    uniqCount=0
    while read begin end; do
      uniqCount=$((uniqCount+end-begin+1))
    done
    uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
    printf '%.5f' "$uniqRatio" > $file.unr
    } < $file
  done

  ## Save result
  echo -e "#mut\tuniqRat\tk" > "result-ancient-A"
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      printf '%d\t%.5f\t%d\n' "$j" "$(< ${i}-A.fa-A_mut_${j}.pos.unr)" "$i" >> "result-ancient-A"
    done
  done
  echo "Finished."

  rm -f *.unr *.sng

  ### Sidron
  echo "Running FRUIT on mutated ancient dataset Sidron ..."

  target=""
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:S_mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r S.fa -t $target -k $i

    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      mv S.fa-S_mut_${j}.sng  ${i}-S.fa-S_mut_${j}.sng;
    done
  done
  
  ## Filter
  target=""
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      target="$target:${i}-S.fa-S_mut_${j}.sng"
    done
  done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-filter -i $target -w 5
  done

  # Remove empty files
  for file in *S.fa*.pos; do
    numLines=$(wc -l < $file)
    if [[ $numLines -eq 1 ]]; then rm -f $file; fi
  done

  echo "Calculating uniqueness ratio ..."
  uniqRatio=0
  for file in *S.fa*.pos; do
    {
    read watermark name size
    uniqCount=0
    while read begin end; do
      uniqCount=$((uniqCount+end-begin+1))
    done
    uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
    printf '%.5f' "$uniqRatio" > $file.unr
    } < $file
  done

  ## Save result
  echo -e "#mut\tuniqRat\tk" > "result-ancient-S"
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      printf '%d\t%.5f\t%d\n' "$j" "$(< ${i}-S.fa-S_mut_${j}.pos.unr)" "$i" >> "result-ancient-S"
    done
  done
  echo "Finished."

  rm -f *.unr *.sng

  ### Vindija
  echo "Running FRUIT on mutated ancient dataset Vindija ..."

  target=""
  for (( i=MIN_MUTE; i<=MAX_MUTE; i++ )); do target="$target:V_mut_$i"; done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-map -r V.fa -t $target -k $i

    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      mv V.fa-V_mut_${j}.sng  ${i}-V.fa-V_mut_${j}.sng;
    done
  done
  
  ## Filter
  target=""
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      target="$target:${i}-V.fa-V_mut_${j}.sng"
    done
  done
  target=${target:1}

  for (( i=MIN_K; i<=MAX_K; i++ )); do
    ./fruit-filter -i $target -w 5
  done

  # Remove empty files
  for file in *V.fa*.pos; do
    numLines=$(wc -l < $file)
    if [[ $numLines -eq 1 ]]; then rm -f $file; fi
  done

  echo "Calculating uniqueness ratio ..."
  uniqRatio=0
  for file in *V.fa*.pos; do
    {
    read watermark name size
    uniqCount=0
    while read begin end; do
      uniqCount=$((uniqCount+end-begin+1))
    done
    uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
    printf '%.5f' "$uniqRatio" > $file.unr
    } < $file
  done

  ## Save result
  echo -e "#mut\tuniqRat\tk" > "result-ancient-V"
  for (( i=MIN_K; i<=MAX_K; i++ )); do
    for (( j=MIN_MUTE; j<=MAX_MUTE; j++ )); do
      printf '%d\t%.5f\t%d\n' "$j" "$(< ${i}-V.fa-V_mut_${j}.pos.unr)" "$i" >> "result-ancient-V"
    done
  done
  echo "Finished."

  rm -f *.unr *.sng

  if [[ ! -d ../../result ]]; then mkdir -p ../../result; fi
  cp result-ancient-A  ../../result
  cp result-ancient-S  ../../result
  cp result-ancient-V  ../../result

  rm -f fruit-map fruit-filter fruit-visual
  cd ../..
fi

########## RUN HS ##########
if [[ $RUN_HS -eq 1 ]]; then
  cp fruit-map fruit-filter fruit-visual dataset/ancient/A.fa \
  dataset/ancient/S.fa dataset/ancient/V.fa  dataset/hs
  cd dataset/hs

  MIN=1
  MAX=20412
  STEP=10000
  TOP=15

  echo "Running FRUIT on Neanderthal and modern human ..."
  ## Map
  for (( i=MIN; i<=MAX; i+=STEP )); do
    target=""
    for (( j=i; j<=MAX && j<i+STEP; j++ )); do target="$target:out$j.fa"; done
    target=${target:1}

    ./fruit-map -r A.fa:S.fa:V.fa -t $target -k 8
  done

  for (( i=MIN; i<=MAX; i++ )); do
    mv out${i}.fa.sng  out${i}.sng;
    rm -f A.fa-out${i}.fa.sng  S.fa-out${i}.fa.sng  V.fa-out${i}.fa.sng
  done

  ### Filter
  for (( i=MIN; i<=MAX; i+=STEP )); do
    target=""
    for (( j=i; j<=MAX && j<i+STEP; j++ )); do
      target="$target:out$j.sng";
    done
    target=${target:1}

    ./fruit-filter -i $target -w 5
  done

  ## Remove empty files
  for file in *.pos; do
    numLines=$(wc -l < $file)
    if [[ $numLines -eq 1 ]]; then rm -f $file; fi
  done

  echo "Calculating uniqueness ratio ..."
  uniqRatio=0
  for file in *.pos; do
    {
    read watermark name size
    uniqCount=0
    while read begin end; do
      uniqCount=$((uniqCount+end-begin+1))
    done
    uniqRatio=$(bc <<< "scale=5; $uniqCount/$size")
    printf '%.5f' "$uniqRatio" > $file.unr
    } < $file
  done

  rm -f uniq
  for file in *.unr; do
    printf '%s\t%.5f\n' "${file/.pos.unr/}" "$(< $file)" >> uniq;
  done
  
  ## Sort uniqueness ratios
  sort -k2,2r -k1,1V uniq > uniq_sorted
  head -n $TOP uniq_sorted > uniq_head
  tail -n $TOP uniq_sorted > uniq_tail
  echo "Finished."

  if [[ ! -d ../../result ]]; then mkdir -p ../../result; fi
  cp uniq*  ../../result

  ### Visual
  target=""
  target=$(awk -F"\t" '{printf ":"$1".pos"}' uniq_head);
  target=${target:1}
  ./fruit-visual -i $target -o head.svg
  
  target=""
  target=$(awk -F"\t" '{printf ":"$1".pos"}' uniq_tail);
  target=${target:1}
  ./fruit-visual -i $target -o tail.svg

  if [[ ! -d ../../result ]]; then mkdir -p ../../result; fi
  mv *.svg  ../../result

  rm -f *.unr  *.sng
  rm -f A.fa S.fa V.fa fruit-*
  cd ../..
fi