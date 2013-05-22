#!/bin/bash
cd "$1"
cp tool-data/nullseq_indices.loc.sample ../../tool-data/nullseq_indices.loc
cp tool-data/sample_roc_chen.png ../../tool-data
cp tool-data/classify_output.out ../../test-data
cp tool-data/classify_test.fa ../../test-data
cp tool-data/kmersvm_output_weights.out ../../test-data
cp tool-data/test_positive.fa ../../test-data
cp tool-data/test_negative.fa ../../test-data
cp tool-data/test_weights.out ../../test-data
cp tool-data/train_predictions.out ../../test-data

