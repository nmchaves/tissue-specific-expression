"""
    Gene Ontology Evidence Codes

"""


class EvidenceCodes:
    def __init__(self):
        self.codes = {}
        # Experimental Evidence
        self.codes['exp'] = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
        # Computational Analysis
        self.codes['compan'] = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA']
        # Author Statement
        self.codes['auth'] = ['TAS', 'NAS']
        # Curatorial
        self.codes['cur'] = ['IC', 'ND']
        # Electronic (Do NOT use for experiments because these are unreliable!)
        self.codes['elec'] = ['IEA']

    def get_codes(self, categories):
        codes = []
        keys = self.codes.keys()
        for category in categories:
            if category in keys:
                codes.extend(self.codes[category])
            else:
                raise category + ' is not a valid category. Use 1 or more of exp, compan, auth, or cur.'
        return codes
