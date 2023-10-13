# TODO
- [ ] Pangolin
- [ ] SpliceAI
- [ ] splam
- [ ] absplice

# Usage
```python
from splicing_models import Pangolin

# sequence I want to score
sequence = "AAAAAAAAAA"
sequence = "AAAAAAAAAC"

# score sequence
model = Pangolin()
model.score_sequence(sequence)
# this gives an error because the sequence is not long enough

# add context
import pyfastx
genome_sequence = pyfastx.Fasta(genome_reference_file)

model = Pangolin()
model.score_sequence(
    sequence, 
    chromosome="chr9", 
    position_start=123456789, 
    position_end=123456789, # this is an insertion 
    strand="-", 
    genome_sequence=genome_sequence
)

# add 5000 upstream and downstream context


# score the sequence

# compare the scores of two sequences
import pyfastx
genome_sequence = pyfastx.Fasta(genome_reference_file)
chromosome="chr9"
position_start=123456789
position_end=123456789 # this is an insertion 
strand="-"
genome_sequence=genome_sequence

sequence_ref = "A"
sequence_alt = "AACGT"
position_alt = 1

model = Pangolin()
scores_ref = model.score_sequence(sequence_ref, strand, chromosome, position_start, position_end, genome_sequence)
scores_alt = model.score_sequence(sequence_alt, strand, chromosome, position_start, position_end, genome_sequence)
loss, gain = model.compare_scores(scores_alt, scores_ref, position_alt)
```


# Published models predicting alternative splicing

- Sequence inputs
  - [SpliceAI](https://github.com/Illumina/SpliceAI):
    - annotates genetic variants with their predicted effect on splicing
    - https://doi.org/10.1016/j.cell.2018.12.015
    - training code: https://github.com/akashc1/splice_2019
  - [DARTS](https://github.com/Xinglab/DARTS):
    - deep-learning Augmented RNA-seq analysis of Transcript Splicing
    - https://doi.org/10.1038/s41592-019-0351-9
  - [Pangolin](https://github.com/tkzeng/Pangolin)
    - deep-learning method for predicting splice site strengths
    - https://doi.org/10.1186/s13059-022-02664-4
  - [splam](https://github.com/Kuanhao-Chao/splam)
    - deep learning splice site predictor that improves spliced alignments
    - https://doi.org/10.1101/2023.07.27.550754
  - [Deep Splicing Code (DSC)](https://github.com/louadi/DSC)
    - deep learning approach for categorizing the well-studied classes of AS namely alternatively skipped exons, alternative 5’ss, alternative 3’ss, and constitutively spliced exons based only on the sequence of the exon junctions
    - https://doi.org/10.3390/genes10080587
  - [EnsembleSplice](https://github.com/OluwadareLab/EnsembleSplice)
    - Ensemble Deep Learning for Splice Site Prediction
    - https://doi.org/10.1186/s12859-022-04971-w
  - [AttentionSplice](https://github.com/EvilBoom/Attention_Splice)
    - AttentionSplice: An interpretable multi-head self-attention-based hybrid deep learning model in splice site prediction
    - https://doi.org/10.1049/cje.2021.00.221
  - [AbSplice](https://github.com/gagneurlab/absplice)
    - aberrant splicing prediction across human tissues
    - https://www.nature.com/articles/s41588-023-01373-3
- [SpliceBERT](https://github.com/biomed-AI/SpliceBERT)
  - Self-supervised learning on millions of pre-mRNA sequences improves sequence-based RNA splicing prediction
  - https://doi.org/10.1101/2023.01.31.526427
- [borzoi](https://github.com/calico/borzoi)
  - Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation
  - https://doi.org/10.1101/2023.08.30.555582
- [SpliceVault](https://github.com/kidsneuro-lab/SpliceVault)
  - SpliceVault predicts the precise nature of variant-associated mis-splicing
  - https://www.nature.com/articles/s41588-022-01293-8
- [BigRNA]()
  - An RNA foundation model enables discovery of disease mechanisms and candidate therapeutics
  - https://doi.org/10.1101/2023.09.20.558508
- [Interpretable Splicing Model](https://github.com/regev-lab/interpretable-splicing-model)
  - Deciphering RNA splicing logic with interpretable machine learning
  - https://doi.org/10.1073/pnas.2221165120

- Methylation inputs
  - [methylAltSplicing](https://github.com/BauerLab/methylAltSplicing)


- Other
  - [Benchmarking_splice_prediction_tools](https://github.com/cmbi/Benchmarking_splice_prediction_tools)
    - benchmarking tool
    - https://doi.org/10.1002/humu.24212
