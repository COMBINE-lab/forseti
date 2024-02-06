# class ForsetiReads:
#     def __init__(self, read_name, is_multi_mapped, forseti_predictions):
#         self.read_name = read_name
#         self.is_multi_mapped = is_multi_mapped
#         self.predictions = forseti_predictions


class ForsetiReads:
    def __init__(self, read_name, is_multi_mapped, forseti_predictions=None):
        self.read_name = read_name
        self.is_multi_mapped = is_multi_mapped
        # Initialize predictions as an empty list if None is provided
        self.predictions = forseti_predictions if forseti_predictions is not None else []

    def add_prediction(self, prediction):
        """Add a ForsetiPredictions object to the predictions list."""
        self.predictions.append(prediction)

    def get_predictions(self):
        """Return the list of ForsetiPredictions objects."""
        return self.predictions

    # If needed, you can add more methods to remove or update predictions


class ForsetiPredictions:
    def __init__(self, gene_id, orientation, splicing_status, max_prob, unsplice_prob, splice_prob):
        self.gene_id = gene_id
        self.orientation = orientation
        self.splicing_status = splicing_status
        self.max_prob = max_prob
        self.unsplice_prob = unsplice_prob
        self.splice_prob = splice_prob
