import argparse
import copy
import dataclasses
import functools
import itertools
import json
import os
import shutil
from typing import (
    Sequence,
)
import plotly.graph_objects as go
import numpy as np
import numpy.typing as npt
import scipy.ndimage
import scipy.fft
import matplotlib.pyplot as plt


def despike_hampel_2d(z, k, win):
    med = scipy.ndimage.median_filter(z, size=win)
    mad = scipy.ndimage.median_filter(np.abs(z - med), size=win) + 1e-12
    # 正規分布仮定のスケール係数 1.4826
    robust_z = (z - med) / (1.4826 * mad)
    out = np.abs(robust_z) > k
    z2 = z.copy()
    z2[out] = med[out]
    return z2, out


@dataclasses.dataclass
class ChevronPatternConfig:
    hampel_filter_threshold: float
    hampel_filter_window_drive_frequency: int
    hampel_filter_window_time: int
    time_spectrum_freq_band_max: int
    peak_strength_thresholds: Sequence[float]

    def __post_init__(self):
        self._validate()

    def _validate(self):
        if self.hampel_filter_threshold <= 0:
            raise ValueError('hampel_filter_threshold must be positive')

        if self.hampel_filter_window_drive_frequency <= 0:
            raise ValueError('hampel_filter_window_drive_frequency must be positive')

        if self.hampel_filter_window_time <= 0:
            raise ValueError('hampel_filter_window_time must be positive')

        if self.time_spectrum_freq_band_max <= 0:
            raise ValueError('time_spectrum_freq_band_max must be positive')

        if len(self.peak_strength_thresholds) == 0 or any(
            b <= a for a, b in itertools.pairwise(self.peak_strength_thresholds)
        ):
            raise ValueError('peak_strength_thresholds must be strictly increasing')


class ChevronPattern:
    def __init__(
        self,
        xs: Sequence[float],
        ys: Sequence[float],
        zs: Sequence[Sequence[float]],
        config: ChevronPatternConfig,
    ):
        self.xs = np.asarray(xs, dtype=np.float64)
        self.ys = np.asarray(ys, dtype=np.float64)
        self.zs = np.asarray(zs, dtype=np.float64)
        self.config = config

        self._validate_input()

    @functools.cached_property
    def _despike_result(self):
        return despike_hampel_2d(
            self.zs,
            k=self.config.hampel_filter_threshold,
            win=(
                self.config.hampel_filter_window_time,
                self.config.hampel_filter_window_drive_frequency,
            ),
        )

    @functools.cached_property
    def zs_despiked(self):
        return self._despike_result[0]

    @functools.cached_property
    def spike_mask(self):
        return self._despike_result[1]

    @functools.cached_property
    def spectrum(self) -> npt.NDArray[np.complex128]:
        res = scipy.fft.rfft2(self.zs_despiked, axes=(1, 0))
        return res

    @functools.cached_property
    def spectrum_folded(self) -> npt.NDArray[np.float64]:
        magnitude = np.abs(self.spectrum)
        magnitude += np.flip(magnitude, axis=1)
        magnitude = magnitude[:, : magnitude.shape[1] // 2]
        return magnitude

    @functools.cached_property
    def spectrum_time(self):
        return np.mean(
            self.spectrum_folded[:, : self.config.time_spectrum_freq_band_max + 1],
            axis=1,
        )

    @functools.cached_property
    def rabi_cycles(self):
        return int(np.argmax(self.spectrum_time[1:]) + 1)

    @functools.cached_property
    def peak_strength(self):
        return self.spectrum_time[self.rabi_cycles] / np.median(self.spectrum_time[1:])

    @functools.cached_property
    def quality_level(self):
        quality_level = np.searchsorted(
            self.config.peak_strength_thresholds, self.peak_strength, side='left'
        )
        return int(quality_level)

    @functools.cached_property
    def num_spikes(self):
        return int(sum(self.spike_mask.ravel()))

    def plot_fft(self):
        magnitude = np.abs(self.spectrum)
        plt.imshow(np.log(1 + magnitude), cmap='viridis')
        plt.title('rfft2 Magnitude Spectrum')
        plt.show()

    def plot_fft_folded(self):
        magnitude = np.abs(self.spectrum_folded)
        plt.imshow(np.log(1 + magnitude), cmap='viridis')
        plt.title('rfft2 Magnitude Spectrum')
        plt.show()

    def save_fft(self, filename):
        magnitude = np.abs(self.spectrum)
        plt.imshow(np.log(1 + magnitude), cmap='viridis')
        plt.title('rfft2 Magnitude Spectrum')
        plt.savefig(filename)

    def save_fft_folded(self, filename):
        plt.imshow(np.log(1 + self.spectrum_folded), cmap='viridis')
        plt.title('rfft2 Magnitude Spectrum')
        plt.savefig(filename)

    def _validate_input(self):
        if self.zs.ndim != 2:
            raise ValueError(f'zs must be 2D, got {self.zs.ndim}D')
        if self.zs.shape != (len(self.ys), len(self.xs)):
            raise ValueError(
                f'shape mismatch: zs{self.zs.shape} vs (len(ys),len(xs))={(len(self.ys),len(self.xs))}'
            )
        if not np.all(np.isfinite(self.zs)):
            raise ValueError('zs contains NaN/Inf')
        if len(self.xs) < 2 or len(self.ys) < 2:
            raise ValueError('xs/ys too short')
        if np.any(np.diff(self.xs) <= 0):
            raise ValueError('xs must be strictly increasing')
        if np.any(np.diff(self.ys) <= 0):
            raise ValueError('ys must be strictly increasing')


def create_figure(data, zs=None):
    if zs is not None:
        data = copy.deepcopy(data)
        data['data'][0]['z'] = zs.tolist()

    return go.Figure(**data)


def process_data(
    data,
    conf,
    image_dir_base=None,
    plot=False,
    json_output=False,
    debug=False,
):
    try:
        chevron_pattern = ChevronPattern(
            data['data'][0]['x'], data['data'][0]['y'], data['data'][0]['z'], conf
        )
    except Exception as e:
        if json_output:
            result = {
                'rabi_cycles': None,
                'quality_level': None,
                'status': 'ERROR',
                'error': str(e),
            }
            if debug:
                result |= {
                    'num_spikes': None,
                    'peak_strength': None,
                }
            print(json.dumps(result))
            return
        raise

    if image_dir_base:
        image_dir = os.path.join(image_dir_base, str(chevron_pattern.quality_level))
        os.makedirs(image_dir, exist_ok=True)
        qubit_idx = data['layout']['title']['text'][-3:]

        fig1 = create_figure(data)
        fig1_path1 = os.path.join(image_dir_base, f'qubit_{qubit_idx}_1_orig.png')
        fig1_path2 = os.path.join(
            image_dir, f'qubit_{qubit_idx}_{chevron_pattern.rabi_cycles:02}.png'
        )
        fig1.write_image(fig1_path1)
        shutil.copy(fig1_path1, fig1_path2)

        fig2 = create_figure(data, chevron_pattern.spike_mask.astype(int))
        fig2.write_image(
            os.path.join(image_dir_base, f'qubit_{qubit_idx}_2_spike_mask.png')
        )

        chevron_pattern.save_fft(
            os.path.join(image_dir_base, f'qubit_{qubit_idx}_3_fft.png')
        )

        chevron_pattern.save_fft_folded(
            os.path.join(image_dir_base, f'qubit_{qubit_idx}_4_fft_folded.png')
        )

        fig5 = create_figure(data, chevron_pattern.zs_despiked.astype(int))
        fig5.write_image(
            os.path.join(image_dir_base, f'qubit_{qubit_idx}_5_despiked.png')
        )

    if plot:
        chevron_pattern.plot_fft_folded()

    if json_output:
        result = {
            'rabi_cycles': chevron_pattern.rabi_cycles,
            'quality_level': chevron_pattern.quality_level,
            'status': 'OK',
            'error': None,
        }
        if debug:
            result |= {
                'num_spikes': chevron_pattern.num_spikes,
                'peak_strength': chevron_pattern.peak_strength,
            }
        print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input-file', required=True)
    parser.add_argument('-c', '--conf-file', required=True)
    parser.add_argument('--image-dir')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--json', action='store_true')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    with open(args.conf_file) as f:
        conf = ChevronPatternConfig(**json.load(f))

    with open(args.input_file) as f:
        data = json.load(f)

    process_data(
        data,
        conf,
        args.image_dir,
        args.plot,
        args.json,
        args.debug,
    )


if __name__ == '__main__':
    main()
