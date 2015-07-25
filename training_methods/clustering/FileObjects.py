__author__ = 'ddeconti'

class DrugEffects():
    def __init__(self, drug_name, se_count, gene_count, smiles):
        self.drug_name = drug_name
        self.side_effect_count = se_count
        self.gene_effect_count = gene_count
        self.SMILES = smiles
        # mutable constructs
        self.rf_prediction_prob = None


    def set_rf_prediction_prob(self, prob):
        self.rf_prediction_prob = prob


    # Getters

    def get_drug_name(self):
        return self.drug_name

    def get_side_effect_count(self):
        return self.side_effect_count

    def get_gene_effect_count(self):
        return self.gene_effect_count

    def get_SMILES(self):
        return self.SMILES

    def get_rf_prediction_prob(self):
        return self.rf_prediction_prob
