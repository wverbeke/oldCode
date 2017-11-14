#!/bin/bash
mkdir -p ../plots/
cd ../plots/
mkdir -p ttgamma
mkdir -p ewkino
mkdir -p tZq
#set up ewkino plots
cd ewkino
mkdir -p dilepCR
cd dilepCR
#make directories for every run
for dir in all2017 RunA RunB RunC RunD RunE RunF
    do mkdir -p $dir
done
for dir in ./*
    do for subdir in inclusive ee em mm
        do mkdir -p $dir/$subdir
    done
done
cd ../..

#set up ttgamma plots
cd ttgamma
for dir in inclusive ee em mm
    do mkdir -p $dir
done
cd ..

#set up tZq plots 
cd tZq 
for dir in inclusive 0bJets_01Jets 0bJets_2Jets 1bJet_01jets 1bJet_23Jets 1bJet_3Jets 2bJets 
    do mkdir -p $dir
done 
cd ..
