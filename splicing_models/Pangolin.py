from pangolin.model import L, W, AR, Pangolin
from pangolin.pangolin import one_hot_encode


class Pangolin:
    def __init__(self):
        """
        - load models
        - 
        """
        
        # load models
        models = []
        for i in [0,2,4,6]:
            for j in range(1,4):
                model = Pangolin(L, W, AR)
                if torch.cuda.is_available():
                    model.cuda()
                    weights = torch.load(resource_filename(__name__,"models/final.%s.%s.3.v2" % (j, i)))
                else:
                    weights = torch.load(resource_filename(__name__,"models/final.%s.%s.3.v2" % (j, i)), map_location=torch.device('cpu'))
                model.load_state_dict(weights)
                model.eval()
                models.append(model)
                
        self.models = models
    
    
    def add_sequence_context(
        self, 
        sequence, 
        chromosome, 
        position_start, 
        position_end, 
        genome_sequence
    ):
        """
        position: position of the first nucleotide in the 
        genome_sequence
        """
        
        sequence = genome_sequence[chromosome][
            int(position_start) - 5001 : int(position_end) + 4999
        ].seq
        
        return sequence
    
    
    def compute_score(self, sequence, strand):
        """
        Outputs scores for each position and 4 tissue models.
        """
        
        sequence = one_hot_encode(sequence, strand).T
        sequence = torch.from_numpy(np.expand_dims(sequence, axis=0)).float()

        if torch.cuda.is_available():
            sequence = sequence.to(torch.device("cuda"))

        scores = []
        for j in range(4):
            scores_tissue = []
            for model in self.models[3*j:3*j+3]:
                with torch.no_grad():
                    score_model_i = model(sequence)[0][[1,4,7,10][j],:].cpu().numpy()
                    if strand == '-':
                        score_model_i = score_model_i[::-1]
                    scores_tissue.append(score_model_i)
            scores.append(np.mean(scores_tissue, axis=0))

        scores = np.array(scores)
        return scores
    
        
    def score_sequence(
        self, 
        sequence, 
        strand,
        chromosome=None, 
        position_start=None, 
        position_end=None, 
        genome_sequence=None
    ):
        """
        get a score for each position of a sequence.
        Note: sequence must have 5000 nt upstream and downstream
        giving context.
        """
        
        # either you give me a sequence with context or context info
        assert (len(sequence)>10000) | (chromosome is not None)
        
        if chromosome is not None:
            # we add context
            sequence = self.add_sequence_context(sequence, chromosome, position_start, position_end, genome_sequence)
        
        # encode sequence
        sequence = one_hot_encode(sequence, strand).T
        
        # score sequence
        scores = compute_score(sequence, strand)
        
        return scores
    
    
    def compare_scores(
        self, 
        scores_alt, 
        scores_ref,
        position_alt
    ):
        """
        compare the scores for each position in two sequences with respect to reference sequence.
        position_alt: 0-based position of where the alteration has happened in the reference sequence.
        """
        
        len_alt = scores_alt.shape[1]
        len_ref = scores_ref.shape[1]
        seq_len_diff = len_alt - len_ref
        
        if seq_len_diff == 0:
            
            # substitutions
            diff = scores_alt - scores_ref
            loss = diff[np.argmin(diff, axis=0), np.arange(diff.shape[1])]
            gain = diff[np.argmax(diff, axis=0), np.arange(diff.shape[1])]
            
        elif seq_len_diff > 0:
            
            # alt had an insertion: do not consider scores in alt insertion region.
            idx = np.arrange(
                set(np.arange(0, len_alt)) - set(np.arange(position_alt, position_alt+seq_len_diff))
            )
            scores_alt = scores_alt[:,idx]
            
            # compute score differences
            diff = scores_alt - scores_ref
            loss = diff[np.argmin(diff, axis=0), np.arange(diff.shape[1])]
            gain = diff[np.argmax(diff, axis=0), np.arange(diff.shape[1])]
            
            
        elif seq_len_diff < 0:
            
            # alt had a deletion: insert 0s in deletion region.
            zeros_to_insert = np.zeros((scores_alt.shape[0], abs(seq_len_diff)))
            scores_alt = np.concatenate(
                [scores_alt[:,:position_alt], zeros_to_insert, scores_alt[:,position_alt:]],
                axis=1
            )
            
            # compute score differences
            diff = scores_alt - scores_ref
            loss = diff[np.argmin(diff, axis=0), np.arange(diff.shape[1])]
            gain = diff[np.argmax(diff, axis=0), np.arange(diff.shape[1])]
            
            
        return loss, gain