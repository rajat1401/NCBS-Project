import random
from typing import List, Dict, Tuple, Iterable


def unroll(xs: Iterable):
    return [x for x in xs]


class Codon:
    def __init__(self, value: str):
        self.value = value
        assert len(value) == 3

    def __eq__(self, other):
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value
    
    def __hash__(self):
        return hash(self.value)

    def __repr__(self):
        return self.value

    def __str__(self):
        return self.value

    def compare(self, other):
        val = list(other.value)
        self_val = list(self.value)
        dist = 0
        for x,y in zip(val, self_val):
            if x != y:
                dist += 1
        return dist


    def get_neighbours(self) -> List:
        val = list(self.value)
        assert len(val) == 3
        nucleotides = [Nucleotide(base) for base in val]
        ret = []
        for i, base in enumerate(nucleotides):
            adj_nucleotides = base.get_others()
            for adj_nucleotide in adj_nucleotides:
                new_nucleotides = nucleotides
                new_nucleotides[i] = adj_nucleotide
                ret.append(get_codon_from_nucleotides(new_nucleotides))
        return ret


class Nucleotide:
    def __init__(self, base: str):
        self.value = base
        assert self.value in ['a', 'u', 'g', 'c']

    def __eq__(self, other):
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value

    def __hash__(self):
        return hash(self.value)

    def __repr__(self):
        return self.value

    def __str__(self):
        return self.value

    def get_others(self) -> List:
        if self.value == 'a':
            return [Nucleotide('u'), Nucleotide('g'), Nucleotide('c')]
        if self.value == 'u':
            return [Nucleotide('a'), Nucleotide('g'), Nucleotide('c')]
        if self.value == 'g':
            return [Nucleotide('u'), Nucleotide('a'), Nucleotide('c')]
        if self.value == 'c':
            return [Nucleotide('u'), Nucleotide('g'), Nucleotide('a')]


def get_codon_from_nucleotides(nucleotides: List[Nucleotide]):
    values = [nucleotide.value for nucleotide in nucleotides]
    return Codon("".join(values))


class Data:
    def __init__(self):
        self.nucleotides = [Nucleotide('a'), Nucleotide('u'), Nucleotide('g'), Nucleotide('c')]
        self.codons = list(set(get_codon_from_nucleotides([nucl1, nucl2, nucl3])
                              for nucl1 in self.nucleotides
                              for nucl2 in self.nucleotides
                              for nucl3 in self.nucleotides))
        assert len(self.codons) == 64

        self.current_genetic_code = {
            Codon('uuu'): 'F',
            Codon('uuc'): 'F',
            Codon('uua'): 'L',
            Codon('uug'): 'L',
            Codon('cuu'): 'L',
            Codon('cuc'): 'L',
            Codon('cua'): 'L',
            Codon('cug'): 'L',
            Codon('auu'): 'I',
            Codon('auc'): 'I',
            Codon('aua'): 'I',
            Codon('aug'): 'M',
            Codon('guu'): 'V',
            Codon('guc'): 'V',
            Codon('gua'): 'V',
            Codon('gug'): 'V',
            Codon('ucu'): 'S',
            Codon('ucc'): 'S',
            Codon('uca'): 'S',
            Codon('ucg'): 'S',
            Codon('ccu'): 'P',
            Codon('ccc'): 'P',
            Codon('cca'): 'P',
            Codon('ccg'): 'P',
            Codon('acu'): 'T',
            Codon('acc'): 'T',
            Codon('aca'): 'T',
            Codon('acg'): 'T',
            Codon('gcu'): 'A',
            Codon('gcc'): 'A',
            Codon('gca'): 'A',
            Codon('gcg'): 'A',
            Codon('uau'): 'Y',
            Codon('uac'): 'Y',
            Codon('uaa'): '*',
            Codon('uag'): '*',
            Codon('cau'): 'H',
            Codon('cac'): 'H',
            Codon('caa'): 'Q',
            Codon('cag'): 'Q',
            Codon('aau'): 'N',
            Codon('aac'): 'N',
            Codon('aaa'): 'K',
            Codon('aag'): 'K',
            Codon('gau'): 'D',
            Codon('gac'): 'D',
            Codon('gaa'): 'E',
            Codon('gag'): 'E',
            Codon('ugu'): 'C',
            Codon('ugc'): 'C',
            Codon('uga'): '*',
            Codon('ugg'): 'W',
            Codon('cgu'): 'R',
            Codon('cgc'): 'R',
            Codon('cga'): 'R',
            Codon('cgg'): 'R',
            Codon('agu'): 'S',
            Codon('agc'): 'S',
            Codon('aga'): 'R',
            Codon('agg'): 'R',
            Codon('ggu'): 'G',
            Codon('ggc'): 'G',
            Codon('gga'): 'G',
            Codon('ggg'): 'G'
        }

        self.outputs = list(self.current_genetic_code.values())

    def compare_codes(self, code1: Dict[Codon, str], code2: Dict[Codon, str]) ->int:
        diff = 0
        # every codon has to be in the domain of every code
        for codon in self.codons:
            if code1[codon] != code2[codon]:
                diff += 1
                # print(codon)
        return diff

    def check_vocabulary(self) -> bool:
        for codon in self.codons:
            if codon not in self.current_genetic_code:
                print("INCONSISTENT VOCABULARY")
                return False
        print("CONSISTENT VOCABULARY")
        return True

    def sampleDistinctCodons(self, n) -> List[Codon]:
        assert n < len(self.codons)
        codon1 = self.codons[random.randint(0, len(self.codons)-1)]
        sampled = set()
        sampled.add(codon1)
        codon2 = self.codons[random.randint(0, len(self.codons)-1)]
        while len(sampled) < n:
            while codon2 in sampled:
                codon2 = self.codons[random.randint(0, len(self.codons)-1)]
            sampled.add(codon2)
        return list(sampled)

    def validate_code(self, code: Dict[Codon, str]) -> bool:
        return len(unroll(code.values())) == len(self.codons)

    def sample_without_swap(self, num_edits: int) -> Dict[Codon, str]:
        # emulating do while
        if num_edits == 0:
            return self.current_genetic_code.copy()
        new_code2 = self.current_genetic_code.copy()
        sampled_codons = self.sampleDistinctCodons(num_edits)
        for codon in sampled_codons:
            prev_code = new_code2.copy()
            new_code2[codon] = self.outputs[random.randint(0, len(self.outputs)-1)]
            if new_code2[codon] == prev_code[codon]:
                print("repetition")
                new_code2[codon] = self.outputs[random.randint(0, len(self.outputs) - 1)]

            while not self.validate_code(new_code2):
                new_code2 = prev_code.copy()
                new_code2[codon] = self.outputs[random.randint(0, len(self.outputs)-1)]
                if new_code2[codon] == prev_code[codon]:
                    print("repetition")
                    new_code2[codon] = self.outputs[random.randint(0, len(self.outputs) - 1)]
        return new_code2

    def sample_1_swap(self, prev_code: Dict[Codon, str]) -> Dict[Codon, str]:
        aa1, aa2 = self.sampleDistinctCodons(2)
        while prev_code[aa1] == prev_code[aa2]:
            aa1, aa2 = self.sampleDistinctCodons(2)
        new_code2 = self.current_genetic_code.copy()
        new_code2[aa1], new_code2[aa2] = new_code2[aa2], new_code2[aa1]
        return new_code2

    def sample_n_swaps(self, num_swaps: int) -> Dict[Codon, str]:
        codon_indices = random.sample(range(len(self.codons)), 2*num_swaps)
        codons = [self.codons[codon_index] for codon_index in codon_indices]
        pairs_codons = [[codons[i], codons[i+1]] for i in range(0, 2*num_swaps, 2)]
        ret = self.current_genetic_code.copy()
        for pair in pairs_codons:
            ret[pair[0]], ret[pair[1]] = ret[pair[1]], ret[pair[0]]
        if self.compare_codes(ret, self.current_genetic_code) < 2*num_swaps:
            print("DISTANCE IS LESS FOR " + str(num_swaps)+" SWAPS")
        return ret

if __name__ == '__main__':
    data = Data()
    new_code = data.sample_without_swap(3)
    print(data.compare_codes(new_code, data.current_genetic_code))
    swap_code = data.sample_n_swaps(4)
    print(data.compare_codes(swap_code, data.current_genetic_code))
    # print(new_code)
    # data.check_vocabulary()
    