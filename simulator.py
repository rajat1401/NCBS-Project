from generator import *
import numpy as np
from typing import List, Dict, Tuple
import pickle
import time

#PROBABILTY OF MISINTERPRETATION IS THE NORMALIZED EDIT DISTANCE
aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

random.seed(0)

def get_misinterpretations(codon: Codon, c2id: Dict[Codon, int]) -> np.ndarray:
    all_codons = unroll(c2id.keys())
    distances = np.zeros(len(all_codons),) # type is float
    for other_codon in all_codons:
        #todo if other_codon==codon then what? check Koonin paper
        distances[c2id[other_codon]] = codon.compare(other_codon)
    #todo distances[codon] = ?
    return distances # shape is (64,)

# 0, 9, 27, 27
def get_kdistance_codons(codon: Codon, c2id: Dict[Codon, int], k: int) -> List[Codon]:
    ls = []
    for c in c2id:
        if (codon.compare(c) ==k):
            ls.append(c)
    print(len(ls))
    return ls

def codons_to_ids(data: Data) -> Dict[Codon, int]:
    c2id = {}
    i = 0
    for codon in data.codons:
        c2id[codon] = i
        i += 1
    return c2id

def aminoacids_to_ids(data: Data) -> Dict[str, int]:
    aa2id = {
        aas[i] : i for i in range(len(aas))
    }

    return aa2id

#WHY EXPONENTIAL, WHY NOT POWER LAW?
def convert_dist_to_prob(dist: np.ndarray, baserate) -> np.ndarray:
    #todo vary baserates

    normalizing_factor = sum(baserate ** (dist))
    mistrans_prob = baserate **(dist)/normalizing_factor
    total_prob = sum(mistrans_prob)
    assert total_prob > 0.99 and total_prob < 1.01, print("total_prob" + str(total_prob))
    # print(unroll(zip(dist,mistrans_prob)))
    return mistrans_prob # shape is (64,)

def get_all_mistrans_dist(c2id: Dict[Codon, int]) -> np.ndarray:
    probs = np.zeros((len(c2id), len(c2id)))
    for codon in c2id:
        prob = get_misinterpretations(codon, c2id)
        probs[c2id[codon],:] = prob
    return probs #shape is (64,64) prob[x,y] is prob that x was misread as y

def get_all_mistrans_prob(c2id: Dict[Codon, int], baserate: float) -> np.ndarray:
    probs = np.zeros((len(c2id), len(c2id)))
    for codon in c2id:
        prob = convert_dist_to_prob(get_misinterpretations(codon, c2id), baserate)
        probs[c2id[codon],:] = prob
    return probs #shape is (64,64)

def get_mistranslation_penalties(aa2id: Dict[str, int]) -> np.ndarray:
    penalties = np.genfromtxt('Blosum62.txt')
    penalties = np.exp(-penalties)
    for i in range(penalties.shape[0]):
        penalties[i,i]=0.

    num_outputs = len(aa2id)
    #make the matrix symmetric
    # for i in range(num_outputs):
    #     for j in range(i, num_outputs):
    #         penalties[i,j] = penalties[i,i] + penalties[j,j] - 2*penalties[i,j]
    #     penalties[i,i] = 0.

    #shape (21,21)
    return penalties

# Measure of how good or bad a genetic code is
def phi(mistrans_prob: np.ndarray, penalties: np.ndarray,
        c2id: Dict[Codon, int], new_code: Dict[Codon, str], aminoacid2id: Dict[str, int]) -> np.ndarray:
    '''

    :param mistrans_prob: 64,64
    :param penalties: 21,21
    :param c2id:
    :param new_code:
    :param aminoacid2id:
    :return:
    '''
    # todo what if some amino acids don't exist in the new genetic code? Account for them
    ret = np.zeros((len(c2id)))
    for codon1 in c2id:
        codon1_id = c2id[codon1]
        total = 0
        for codon2 in c2id:
            codon2_id = c2id[codon2]
            aa1_id =  aminoacid2id[new_code[codon1]]
            aa2_id =  aminoacid2id[new_code[codon2]]
            update = mistrans_prob[codon1_id, codon2_id] * penalties[aa1_id, aa2_id]
            if codon1_id != codon2_id:
                total += update # todo why this?
        ret[codon1_id] = total
    return ret # shape is (64,)

def gaussian(mean: np.ndarray, variance: np.ndarray) -> np.ndarray:

    covariance = np.eye(mean.shape[0]) * variance
    new_xs = np.random.multivariate_normal(mean, covariance)
    return new_xs

def random_walk_aminoacids(num_walks: int, aa2id: Dict[str, int]):
    proteins = []
    acids = [acid for acid in aa2id.keys() if acid != "*"]
    for walk in range(num_walks):
        protein = []
        protein.append("Met")
        for _ in range(98):
            index = random.randint(0, len(acids)-1)
            protein.append(acids[index])
        protein.append("*")
        proteins.append(protein)
    with open("proteins.pkl","wb") as f:
        pickle.dump(proteins,f)
    return




def normalize_populations(vec: np.ndarray) -> np.ndarray:
    #todo should we do this or use Monod constant or something like that
    vec[vec < 0.] = 0.
    return (vec/np.sum(vec))*100

def get_variance(curr_populations: np.ndarray):
    return curr_populations * np.ceil(10./(curr_populations+10))

def convert_vector_to_str(arr: np.ndarray):
    size = arr.shape[0]
    assert len(arr.shape) == 1
    ret = ""
    for k in range(size):
        ret = ret+str(arr[k])+" "
    ret = ret[:-1]
    ret += "\n"
    return ret

def perform_simulation(init_populations: np.ndarray, k1: float,
             k2: float, k3: float,
             err_rates: np.ndarray, num_generations: int, distances: np.ndarray):
    curr_populations = init_populations
    # print(err_rates)
    # print(curr_populations)
    # print(distances)
    f = open('sample','w')

    f.write(convert_vector_to_str(distances))
    #f.close()
    #time.sleep(2)
    for t in range(num_generations):
        #f.open('sample', 'w')
        delta_x = gaussian(curr_populations * k1, get_variance(curr_populations)) \
                  + k2*(curr_populations*curr_populations)\
                  - k3*err_rates
        # print("************************************************************")
        # print(delta_x)
        curr_populations = normalize_populations(curr_populations + delta_x)
        # plot(curr_populations, distances)
        print(curr_populations)
        f.write(convert_vector_to_str(curr_populations))
        #f.close()
        #time.sleep(1.0)
    f.close()
    print(distances)
    print(np.ceil(err_rates))
    print(init_populations)
    return curr_populations



def prepare_simulation(data: Data, rel_error: float, k1: float, k2: float, k3: float):
    minimum, maximum = 5, 10
    num_codes = random.randint(minimum, maximum)
    #todo change num_codes again
    mean = 100. / num_codes
    init_populations = normalize_populations(mean + np.ceil(np.random.randn(num_codes)))
    translation_mechanism_baserate = {
     1. : 0.08664034996495773,
     10. : 0.40799336352310517,
     0.1 : 0.01076003848545575
    }
    assert isinstance(rel_error, float) and rel_error in translation_mechanism_baserate
    baserate = translation_mechanism_baserate[rel_error]
    range_dist = unroll(range(2,15))
    distances = np.array([random.choice(range_dist) for _ in range(num_codes)])

    #todo change this
    distances[0]=1


    c2id = codons_to_ids(data)
    aa2id = aminoacids_to_ids(data)
    print(distances.shape)
    codes = []
    err_rates = []
    for j, distance in enumerate(distances):
        code = data.sample_without_swap(distance)
        codes.append(code)
        if(distance != data.compare_codes(code, data.current_genetic_code)):
            distances[j] = data.compare_codes(code, data.current_genetic_code)
        penalties = get_mistranslation_penalties(aa2id)
        phi_values = phi(get_all_mistrans_prob(c2id, baserate), penalties, c2id, code, aa2id)
        err_rates.append(np.sum(phi_values))
    # print(err_rates)
    perform_simulation(init_populations, k1, k2, k3, np.array(err_rates), 500, distances)


if __name__ == '__main__':

    k1 = 0.5
    k2 = 0.001
    k3 = 0.01
    k4 = 1.0
    codon = Codon('uua')
    data2 = Data()
    c2id = codons_to_ids(data2)
    aa2id = aminoacids_to_ids(data2)
    get_kdistance_codons(codon, c2id, 0)
    prepare_simulation(data2, 10., k1, k2, k3)
