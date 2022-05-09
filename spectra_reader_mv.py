import subject as sbj
import numpy as np
import os


def info_header(espectro):
    info = {}
    header = all_info_header(espectro)
    header = header.split("\n")
    header = [line.replace('!', '#', 1) for line in header]
    header = [line.replace(':', ' = ', 1) for line in header]
    info['examination_name'] = info_header_formatted(header, 26)
    info['scan_id'] = info_header_formatted(header, 28)
    info['scan_date'] = info_header_formatted(header, 30)
    info['patient_name'] = info_header_formatted(header, 32)
    info['patient_birth_date'] = info_header_formatted(header, 34)
    info['patient_position'] = info_header_formatted(header, 36)
    info['patient_orientation'] = info_header_formatted(header, 38)
    info['samples'] = info_header_formatted(header, 40)
    info['rows'] = info_header_formatted(header, 42)
    info['synthesizer_frequency'] = info_header_formatted(header, 44)
    info['offset_frequency'] = info_header_formatted(header, 46)
    info['sample_frequency'] = info_header_formatted(header, 48)
    info['echo_nr'] = info_header_formatted(header, 50)
    info['nucleus'] = info_header_formatted(header, 54)
    info['t0_mu1_direction'] = info_header_formatted(header, 56)
    info['echo_time'] = info_header_formatted(header, 58)
    info['averages'] = info_header_formatted(header, 62)
    info['volume_selection_enable'] = info_header_formatted(header, 64)
    info['volumes'] = info_header_formatted(header, 66)
    info['ap_size'] = info_header_formatted(header, 68)
    info['lr_size'] = info_header_formatted(header, 70)
    info['cc_size'] = info_header_formatted(header, 72)
    info['ap_off_center'] = info_header_formatted(header, 74)
    info['lr_off_center'] = info_header_formatted(header, 76)
    info['cc_off_center'] = info_header_formatted(header, 78)
    info['ap_angulation'] = info_header_formatted(header, 80)
    info['lr_angulation'] = info_header_formatted(header, 82)
    info['cc_angulation'] = info_header_formatted(header, 84)
    info['volume_selection_method'] = info_header_formatted(header, 86)
    info['t1_measurement_enable'] = info_header_formatted(header, 88)
    info['t2_measurement_enable'] = info_header_formatted(header, 90)
    info['time_series_enable'] = info_header_formatted(header, 92)
    info['phase_encoding_enable'] = info_header_formatted(header, 94)
    info['nr_phase_encoding_profiles'] = info_header_formatted(header, 96)
    info['si_ap_off_center'] = info_header_formatted(header, 98)
    info['si_lr_off_center'] = info_header_formatted(header, 100)
    info['si_cc_off_center'] = info_header_formatted(header, 102)
    info['si_ap_off_angulation'] = info_header_formatted(header, 104)
    info['si_lr_off_angulation'] = info_header_formatted(header, 106)
    info['si_cc_off_angulation'] = info_header_formatted(header, 108)
    info['t0_kx_direction'] = info_header_formatted(header, 110)
    info['t0_ky_direction'] = info_header_formatted(header, 112)
    info['nr_of_phase_encoding_profiles_ky'] = info_header_formatted(header, 114)
    info['phase_encoding_direction'] = info_header_formatted(header, 116)
    info['phase_encoding_fov'] = info_header_formatted(header, 118)
    info['slice_thickness'] = info_header_formatted(header, 120)
    info['image_plane_slice_thickness'] = info_header_formatted(header, 122)
    info['slice_distance'] = info_header_formatted(header, 124)
    info['nr_of_slices_for_multislice'] = info_header_formatted(header, 126)
    info['spec_image_in_plane_transf'] = info_header_formatted(header, 128)
    info['spec_num_col'] = info_header_formatted(header, 142)
    info['spec_col_lower_val'] = info_header_formatted(header, 144)
    info['spec_col_upper_val'] = info_header_formatted(header, 146)
    info['spec_num_row'] = info_header_formatted(header, 156)
    info['spec_row_lower_val'] = info_header_formatted(header, 158)
    info['spec_row_upper_val'] = info_header_formatted(header, 160)
    info['num_dimensions'] = info_header_formatted(header, 172)
    info['dim2_pnts'] = info_header_formatted(header, 192)
    info['dim3_pnts'] = info_header_formatted(header, 206)
    return info


def all_info_header(path_espectro):
    header = open(path_espectro, "r").read()
    return header


def info_header_formatted(header, index):
    if len(header[index].split('=')) == 2:
        return header[index].split('=')[1].replace('\r', '').strip()
    else:
        return ''


# Funcoes auxiliares
# Converts a float in Vax format to IEEE format.
def vax_to_ieee_single_float(data):
    f = []
    nfloat = int(len(data) / 4)
    for i in range(nfloat):
        byte2 = data[0 + i * 4]
        byte1 = data[1 + i * 4]
        byte4 = data[2 + i * 4]
        byte3 = data[3 + i * 4]
        sign = (ord(chr(byte1)) & 0x80) >> 7
        expon = ((ord(chr(byte1)) & 0x7f) << 1) + ((ord(chr(byte2)) & 0x80) >> 7)
        fract = ((ord(chr(byte2)) & 0x7f) << 16) + (ord(chr(byte3)) << 8) + ord(chr(byte4))
        if sign == 0:
            sign_mult = 1.0
        else:
            sign_mult = -1.0;
        if 0 < expon:
            val = sign_mult * (0.5 + (fract / 16777216.0)) * pow(2.0, expon - 128.0)
            f.append(val)
        elif expon == 0 and sign == 0:
            f.append(0)
        else:
            f.append(0)
    return f


def collapse_complexes(data):
    data_iter = iter(data)
    return [complex(r, i) for r, i in zip(data_iter, data_iter)]


# Given a string of data in raw format, returns an iterable (tuple or list) of Python complexes.
def decode_raw(data):
    data = vax_to_ieee_single_float(data)
    data = collapse_complexes(data)
    return data


# data = open('data/s2D_PRESS_144_CC_SENSE_11_1_raw_act.SPAR', "rb").read()


def reader_philips_spar(filename):
    data = open(filename, "rb").read()
    # data = open(filename, "r").read()
    data = decode_raw(data)
    # data is complex64 per Philips documentation for SDAT
    data = np.fromiter(data, np.complex64)
    data = np.conjugate(data)
    return data


# le os arquivos da pasta do param. e monta o objeto subject
def read_filenames(pasta, indiv):

    # lista de individuos
    indivs = [nome for nome in os.listdir(pasta)]
    # lista de dos caminhos das pastas dos individuos
    paths = [os.path.join(pasta, nome) for nome in indivs]

    # monta todos os caminhos dos arquivos de input para o individuo passado como parametro
    arqvs_inputs = [nome for nome in os.listdir(paths[indivs.index(indiv)] + '/inputs')]
    paths_inputs = [os.path.join(paths[indivs.index(indiv)] + '/inputs', nome) for nome in arqvs_inputs]

    arqvs_qm = [nome for nome in os.listdir(paths[indivs.index(indiv)] + '/quality_metrics')]
    paths_qm = [os.path.join(paths[indivs.index(indiv)] + '/quality_metrics', nome) for nome in arqvs_qm]

    # encontra cada arquivo separadamente para retornar na funcao
    spar = [arq for arq in paths_inputs if arq.lower().endswith(".spar")]
    sdat = [arq for arq in paths_inputs if arq.lower().endswith(".sdat")]
    snr = [arq for arq in paths_qm if arq.lower().endswith("snr.txt")]
    fwhmhz = [arq for arq in paths_qm if arq.lower().endswith("fwhmhz.txt")]
    qntf = paths[indivs.index(indiv)] + '/inputs/' + indiv + '.txt'

    return sbj.Subject(spar[0], sdat[0], qntf, snr[0], fwhmhz[0])
