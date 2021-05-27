import math
import bisect
from collections import namedtuple

Match = namedtuple('Match', ['peak1', 'peak2', 'score'])
Alignment = namedtuple('Alignment', ['peak1', 'peak2'])
Peak = namedtuple('Peak',['mz','intensity'])

def normalize_spectrum(spectrum):
    output_spectrum = []
    intermediate_output_spectrum = []
    acc_norm = 0.0
    for s in spectrum:
        intermediate_output_spectrum.append(Peak(s.mz,s.intensity))
        acc_norm += s.intensity**2
    normed_value = math.sqrt(acc_norm)
    for s in intermediate_output_spectrum:
        output_spectrum.append(Peak(s.mz,s.intensity/normed_value))
    return output_spectrum

def sqrt_normalize_spectrum(spectrum):
    output_spectrum = []
    intermediate_output_spectrum = []
    acc_norm = 0.0
    for s in spectrum:
        intensity = math.sqrt(s.intensity)
        intermediate_output_spectrum.append(Peak(s.mz,intensity))
        acc_norm += s.intensity
    normed_value = math.sqrt(acc_norm)
    for s in intermediate_output_spectrum:
        output_spectrum.append(Peak(s.mz,s.intensity/normed_value))
    return output_spectrum

def find_match_peaks_efficient(spec1, spec2, shift, tolerance):
    try:
        return find_match_peaks_mz(spec1, spec2, shift, tolerance)
    except:
        return find_match_peaks_ions(spec1, spec2)

#reimplementation of find_match_peaks, but much more efficient
# Assumes that shift is equal to spec1 - spec2
def find_match_peaks_mz(spec1, spec2, shift, tolerance):
    adj_tolerance =  tolerance + 0.000001
    spec2_mass_list = []

    for i,peak in enumerate(spec2):
        spec2_mass_list.append(peak.mz)

    alignment_mapping = []

    for i, peak in enumerate(spec1):
        left_mz_bound = peak.mz - shift - adj_tolerance
        right_mz_bound = peak.mz - shift + adj_tolerance

        left_bound_index = bisect.bisect_left(spec2_mass_list, left_mz_bound)
        right_bound_index = bisect.bisect_right(spec2_mass_list, right_mz_bound)

        for j in range(left_bound_index,right_bound_index):
            alignment_mapping.append(Alignment(i,j))

    return alignment_mapping

# Assumes that shift is equal to spec1 - spec2
def find_match_peaks_ions(spec1, spec2):

    alignment_mapping = []

    spec1_dict = {ion:i for i,(ion,intensity) in enumerate(spec1)}
    spec2_dict = {ion:i for i,(ion,intensity) in enumerate(spec2)}

    for ion, i in spec1_dict.items():
        if ion in spec2_dict:
            j = spec2_dict[ion]
            alignment_mapping.append(Alignment(i,j))

    return alignment_mapping


# Assumes that shift is equal to spec1 - spec2
def find_match_peaks(spec1,spec2,shift,tolerance):
    adj_tolerance =  tolerance + 0.000001
    low = 0
    high = 0
    alignment_mapping = []
    for i,s1 in enumerate(spec1):
        low = len(spec2)-1
        while low > 0 and (s1.mz - adj_tolerance) < (spec2[low].mz + shift):
            low = low - 1
        while low < len(spec2) and (s1.mz - adj_tolerance) > (spec2[low].mz + shift):
            low = low + 1
        while high < len(spec2) and (s1.mz + adj_tolerance) >= (spec2[high].mz + shift):
            high = high + 1
        for j in range(low,high):
            alignment_mapping.append(Alignment(i,j))
    return alignment_mapping

def alignment_to_match(spec1_n,spec2_n,alignment):
    s1_peak = spec1_n[alignment.peak1].intensity
    s2_peak = spec2_n[alignment.peak2].intensity
    match_score = s1_peak * s2_peak
    return Match(
            peak1 = alignment.peak1,
            peak2 = alignment.peak2,
            score = match_score
    )

###
# This score alignment code will take two spectra
# then it will align the two spectrum given their parent masses
# These spectra are expected to be a list of lists (size two mass, intensity) or a list of tuples
###
def score_alignment(spec1_n, spec2_n, pm1, pm2, tolerance, max_charge_consideration=1, normalize=False, sqrt=False):
    if len(spec1_n) == 0 or len(spec2_n) == 0:
        return 0.0, []

    if normalize:
        if sqrt:
            spec1_n = normalize_spectrum(convert_to_peaks(spec1))
            spec2_n = normalize_spectrum(convert_to_peaks(spec2))
        else:
            spec1_n = normalize_spectrum(convert_to_peaks(spec1))
            spec2_n = normalize_spectrum(convert_to_peaks(spec2))

    shift = (pm1 - pm2)

    zero_shift_alignments = find_match_peaks_efficient(spec1_n,spec2_n,0,tolerance)
    real_shift_alignments = []
    if abs(shift) > tolerance:
        real_shift_alignments = find_match_peaks_efficient(spec1_n,spec2_n,shift,tolerance)

        if max_charge_consideration > 1:
            for charge_considered in range(2, max_charge_consideration + 1):
                real_shift_alignments += find_match_peaks_efficient(spec1_n,spec2_n,shift/charge_considered,tolerance)

    #Making real_shift_alignments without repetition
    real_shift_alignments = list(set(real_shift_alignments))

    zero_shift_match = [alignment_to_match(spec1_n,spec2_n,alignment) for alignment in zero_shift_alignments]
    real_shift_match = [alignment_to_match(spec1_n,spec2_n,alignment) for alignment in real_shift_alignments]

    all_possible_match_scores = zero_shift_match + real_shift_match
    all_possible_match_scores.sort(key=lambda x: x.score, reverse=True)

    reported_alignments = []

    spec1_peak_used = set()
    spec2_peak_used = set()

    total_score = 0.0

    for match in all_possible_match_scores:
        if not match.peak1 in spec1_peak_used and not match.peak2 in spec2_peak_used:
            spec1_peak_used.add(match.peak1)
            spec2_peak_used.add(match.peak2)
            reported_alignments.append(Alignment(match.peak1,match.peak2))
            total_score += match.score

    return total_score, reported_alignments

def score_alignment_matched_peaks(spec1,spec2,pm1,pm2,tolerance,max_charge_consideration=1, reported_alignments=None):
    if reported_alignments == None:
        total_score, reported_alignments = score_alignment_matched_peaks(spec1,spec2,pm1,pm2,tolerance,max_charge_consideration)

    spec1_peak_list = []
    spec2_peak_list = []

    for reported_alignment in reported_alignments:
        spec1_peak_list.append(spec1[reported_alignment.peak1])
        spec2_peak_list.append(spec2[reported_alignment.peak2])

    spec1_n = sqrt_normalize_spectrum(convert_to_peaks(spec1_peak_list))
    spec2_n = sqrt_normalize_spectrum(convert_to_peaks(spec2_peak_list))

    score_total = 0.0
    for i in range(len(spec1_peak_list)):
        score_total += spec1_n[i].intensity * spec2_n[i].intensity

    return score_total
