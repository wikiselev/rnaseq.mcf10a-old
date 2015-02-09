# P110alpha knock in

mv ../H1047R_H-3-EGF/lane2_GTGGCC_L002_R1.fastq.gz ../H1047R_H-3-EGF/ki_3_0.fastq.gz
mv ../H1047R_H-3-EGF/lane2_GTTTCG_L002_R1.fastq.gz ../H1047R_H-3-EGF/ki_3_15.fastq.gz
mv ../H1047R_H-3-EGF/lane2_CGTACG_L002_R1.fastq.gz ../H1047R_H-3-EGF/ki_3_40.fastq.gz
mv ../H1047R_H-3-EGF/lane2_GAGTGG_L002_R1.fastq.gz ../H1047R_H-3-EGF/ki_3_90.fastq.gz
mv ../H1047R_H-3-EGF/lane2_ACTGAT_L002_R1.fastq.gz ../H1047R_H-3-EGF/ki_3_180.fastq.gz
mv ../H1047R_H-3-EGF/lane2_ATTCCT_L002_R1.fastq.gz ../H1047R_H-3-EGF/ki_3_300.fastq.gz

mv ../H1047R_H-2-EGF/lane8_GTGGCC_L008_R1.fastq.gz ../H1047R_H-2-EGF/ki_2_0.fastq.gz
mv ../H1047R_H-2-EGF/lane8_GTTTCG_L008_R1.fastq.gz ../H1047R_H-2-EGF/ki_2_15.fastq.gz
mv ../H1047R_H-2-EGF/lane8_CGTACG_L008_R1.fastq.gz ../H1047R_H-2-EGF/ki_2_40.fastq.gz
mv ../H1047R_H-2-EGF/lane8_GAGTGG_L008_R1.fastq.gz ../H1047R_H-2-EGF/ki_2_90.fastq.gz
mv ../H1047R_H-2-EGF/lane8_ACTGAT_L008_R1.fastq.gz ../H1047R_H-2-EGF/ki_2_180.fastq.gz
mv ../H1047R_H-2-EGF/lane8_ATTCCT_L008_R1.fastq.gz ../H1047R_H-2-EGF/ki_2_300.fastq.gz

mv ../H1047R_H-1-EGF/lane7_GTGGCC_L007_R1.fastq.gz ../H1047R_H-1-EGF/ki_1_0.fastq.gz
mv ../H1047R_H-1-EGF/lane7_GTTTCG_L007_R1.fastq.gz ../H1047R_H-1-EGF/ki_1_15.fastq.gz
mv ../H1047R_H-1-EGF/lane7_CGTACG_L007_R1.fastq.gz ../H1047R_H-1-EGF/ki_1_40.fastq.gz
mv ../H1047R_H-1-EGF/lane7_GAGTGG_L007_R1.fastq.gz ../H1047R_H-1-EGF/ki_1_90.fastq.gz
mv ../H1047R_H-1-EGF/lane7_ACTGAT_L007_R1.fastq.gz ../H1047R_H-1-EGF/ki_1_180.fastq.gz
mv ../H1047R_H-1-EGF/lane7_ATTCCT_L007_R1.fastq.gz ../H1047R_H-1-EGF/ki_1_300.fastq.gz

# PTEN knock out

mv ../PTEN_A-4-EGF/lane4_ATCACG_L004_R1.fastq.gz ../PTEN_A-4-EGF/pten_4_0.fastq.gz
mv ../PTEN_A-4-EGF/lane4_TTAGGC_L004_R1.fastq.gz ../PTEN_A-4-EGF/pten_4_15.fastq.gz
mv ../PTEN_A-4-EGF/lane4_ACTTGA_L004_R1.fastq.gz ../PTEN_A-4-EGF/pten_4_40.fastq.gz
mv ../PTEN_A-4-EGF/lane4_GATCAG_L004_R1.fastq.gz ../PTEN_A-4-EGF/pten_4_90.fastq.gz
mv ../PTEN_A-4-EGF/lane4_TAGCTT_L004_R1.fastq.gz ../PTEN_A-4-EGF/pten_4_180.fastq.gz
mv ../PTEN_A-4-EGF/lane4_GGCTAC_L004_R1.fastq.gz ../PTEN_A-4-EGF/pten_4_300.fastq.gz

mv ../PTEN_A-3-EGF/lane3_ATCACG_L003_R1.fastq.gz ../PTEN_A-3-EGF/pten_3_0.fastq.gz
mv ../PTEN_A-3-EGF/lane3_TTAGGC_L003_R1.fastq.gz ../PTEN_A-3-EGF/pten_3_15.fastq.gz
mv ../PTEN_A-3-EGF/lane3_ACTTGA_L003_R1.fastq.gz ../PTEN_A-3-EGF/pten_3_40.fastq.gz
mv ../PTEN_A-3-EGF/lane3_GATCAG_L003_R1.fastq.gz ../PTEN_A-3-EGF/pten_3_90.fastq.gz
mv ../PTEN_A-3-EGF/lane3_TAGCTT_L003_R1.fastq.gz ../PTEN_A-3-EGF/pten_3_180.fastq.gz
mv ../PTEN_A-3-EGF/lane3_GGCTAC_L003_R1.fastq.gz ../PTEN_A-3-EGF/pten_3_300.fastq.gz

# split the total raw data file into barcode files
perl generic_barcode_splitter_SE ../PTEN_A-2-EGF/lane8_NoIndex_L008_R1.fastq.gz ../PTEN_A-2-EGF/lane8_NoIndex_L008_R2.fastq.gz

mv ../PTEN_A-2-EGF/lane8_ATCACG_L008_R1.fastq.gz ../PTEN_A-2-EGF/pten_2_0.fastq.gz
mv ../PTEN_A-2-EGF/lane8_TTAGGC_L008_R1.fastq.gz ../PTEN_A-2-EGF/pten_2_15.fastq.gz
mv ../PTEN_A-2-EGF/lane8_ACTTGA_L008_R1.fastq.gz ../PTEN_A-2-EGF/pten_2_40.fastq.gz
mv ../PTEN_A-2-EGF/lane8_GATCAG_L008_R1.fastq.gz ../PTEN_A-2-EGF/pten_2_90.fastq.gz
mv ../PTEN_A-2-EGF/lane8_TAGCTT_L008_R1.fastq.gz ../PTEN_A-2-EGF/pten_2_180.fastq.gz
mv ../PTEN_A-2-EGF/lane8_GGCTAC_L008_R1.fastq.gz ../PTEN_A-2-EGF/pten_2_300.fastq.gz

# P110alpha knock out

# split the total raw data file into barcode files
perl generic_barcode_splitter_SE ../INK-4-EGF/lane7_NoIndex_L007_R1.fastq.gz ../INK-4-EGF/lane7_NoIndex_L007_R2.fastq.gz

mv ../INK-4-EGF/lane7_AGTCAA_L007_R1.fastq.gz ../INK-4-EGF/ko_4_0.fastq.gz
mv ../INK-4-EGF/lane7_AGTTCC_L007_R1.fastq.gz ../INK-4-EGF/ko_4_15.fastq.gz
mv ../INK-4-EGF/lane7_ATGTCA_L007_R1.fastq.gz ../INK-4-EGF/ko_4_40.fastq.gz
mv ../INK-4-EGF/lane7_CCGTCC_L007_R1.fastq.gz ../INK-4-EGF/ko_4_90.fastq.gz
mv ../INK-4-EGF/lane7_GTCCGC_L007_R1.fastq.gz ../INK-4-EGF/ko_4_180.fastq.gz
mv ../INK-4-EGF/lane7_GTGAAA_L007_R1.fastq.gz ../INK-4-EGF/ko_4_300.fastq.gz

perl generic_barcode_splitter_SE ../INK-3-EGF/lane6_NoIndex_L006_R1.fastq.gz ../INK-3-EGF/lane6_NoIndex_L006_R2.fastq.gz

mv ../INK-3-EGF/lane6_AGTCAA_L006_R1.fastq.gz ../INK-3-EGF/ko_3_0.fastq.gz
mv ../INK-3-EGF/lane6_AGTTCC_L006_R1.fastq.gz ../INK-3-EGF/ko_3_15.fastq.gz
mv ../INK-3-EGF/lane6_ATGTCA_L006_R1.fastq.gz ../INK-3-EGF/ko_3_40.fastq.gz
mv ../INK-3-EGF/lane6_CCGTCC_L006_R1.fastq.gz ../INK-3-EGF/ko_3_90.fastq.gz
mv ../INK-3-EGF/lane6_GTCCGC_L006_R1.fastq.gz ../INK-3-EGF/ko_3_180.fastq.gz
mv ../INK-3-EGF/lane6_GTGAAA_L006_R1.fastq.gz ../INK-3-EGF/ko_3_300.fastq.gz

mv ../INK-2-EGF/lane2_AGTCAA_L002_R1.fastq.gz ../INK-2-EGF/ko_2_0.fastq.gz
mv ../INK-2-EGF/lane2_AGTTCC_L002_R1.fastq.gz ../INK-2-EGF/ko_2_15.fastq.gz
mv ../INK-2-EGF/lane2_ATGTCA_L002_R1.fastq.gz ../INK-2-EGF/ko_2_40.fastq.gz
mv ../INK-2-EGF/lane2_CCGTCC_L002_R1.fastq.gz ../INK-2-EGF/ko_2_90.fastq.gz
mv ../INK-2-EGF/lane2_GTCCGC_L002_R1.fastq.gz ../INK-2-EGF/ko_2_180.fastq.gz
mv ../INK-2-EGF/lane2_GTGAAA_L002_R1.fastq.gz ../INK-2-EGF/ko_2_300.fastq.gz

# WT

mv ../WT-4-EGF/lane1_CGATGT_L001_R1.fastq.gz ../WT-4-EGF/wt_4_0.fastq.gz
mv ../WT-4-EGF/lane1_TGACCA_L001_R1.fastq.gz ../WT-4-EGF/wt_4_15.fastq.gz
mv ../WT-4-EGF/lane1_ACAGTG_L001_R1.fastq.gz ../WT-4-EGF/wt_4_40.fastq.gz
mv ../WT-4-EGF/lane1_GCCAAT_L001_R1.fastq.gz ../WT-4-EGF/wt_4_90.fastq.gz
mv ../WT-4-EGF/lane1_CAGATC_L001_R1.fastq.gz ../WT-4-EGF/wt_4_180.fastq.gz
mv ../WT-4-EGF/lane1_CTTGTA_L001_R1.fastq.gz ../WT-4-EGF/wt_4_300.fastq.gz

perl generic_barcode_splitter_SE ../WT-3-EGF/lane5_NoIndex_L005_R1.fastq.gz ../WT-3-EGF/lane5_NoIndex_L005_R2.fastq.gz

mv ../WT-3-EGF/lane5_CGATGT_L005_R1.fastq.gz ../WT-3-EGF/wt_3_0.fastq.gz
mv ../WT-3-EGF/lane5_TGACCA_L005_R1.fastq.gz ../WT-3-EGF/wt_3_15.fastq.gz
mv ../WT-3-EGF/lane5_ACAGTG_L005_R1.fastq.gz ../WT-3-EGF/wt_3_40.fastq.gz
mv ../WT-3-EGF/lane5_GCCAAT_L005_R1.fastq.gz ../WT-3-EGF/wt_3_90.fastq.gz
mv ../WT-3-EGF/lane5_CAGATC_L005_R1.fastq.gz ../WT-3-EGF/wt_3_180.fastq.gz
mv ../WT-3-EGF/lane5_CTTGTA_L005_R1.fastq.gz ../WT-3-EGF/wt_3_300.fastq.gz

perl generic_barcode_splitter_SE ../WT-2-EGF/lane4_NoIndex_L004_R1.fastq.gz ../WT-2-EGF/lane4_NoIndex_L004_R2.fastq.gz

mv ../WT-2-EGF/lane4_CGATGT_L004_R1.fastq.gz ../WT-2-EGF/wt_2_0.fastq.gz
mv ../WT-2-EGF/lane4_TGACCA_L004_R1.fastq.gz ../WT-2-EGF/wt_2_15.fastq.gz
mv ../WT-2-EGF/lane4_ACAGTG_L004_R1.fastq.gz ../WT-2-EGF/wt_2_40.fastq.gz
mv ../WT-2-EGF/lane4_GCCAAT_L004_R1.fastq.gz ../WT-2-EGF/wt_2_90.fastq.gz
mv ../WT-2-EGF/lane4_CAGATC_L004_R1.fastq.gz ../WT-2-EGF/wt_2_180.fastq.gz
mv ../WT-2-EGF/lane4_CTTGTA_L004_R1.fastq.gz ../WT-2-EGF/wt_2_300.fastq.gz

# P110alpha knock out + no stimulation

mv ../CT_neg_INK300/lane3_ATCACG_L003_R1.fastq.gz ../CT_neg_INK300/konost_1_300.fastq.gz
mv ../CT_neg_INK300/lane3_TAGCTT_L003_R1.fastq.gz ../CT_neg_INK300/konost_2_300.fastq.gz
mv ../CT_neg_INK300/lane3_GTGGCC_L003_R1.fastq.gz ../CT_neg_INK300/konost_3_300.fastq.gz
mv ../CT_neg_INK300/lane3_TTAGGC_L003_R1.fastq.gz ../CT_neg_INK300/wtb_3_0.fastq.gz

