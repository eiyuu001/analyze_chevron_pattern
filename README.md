## 処理の流れ

1. Hampelフィルタを用いてスパイクノイズを除去. ウィンドウサイズ: `[hampel_filter_window_time, hampel_filter_window_drive_frequency]`, 閾値: `hampel_filter_threshold`
2. 二次元フーリエ変換によりスペクトル情報に変換.
3. スペクトルを入力周波数方向に0 ~ `time_spectrum_freq_band_max`(この値を含む)の範囲で平均し, 時間方向のスペクトルを計算.
4. 3からDC成分(波数0成分)を除去.
5. 4のピーク位置をChevron Patternの時間方向の波数(`rabi_cycles`)とみなす.
6. 4のピーク位置の成分の, 4の中央値に対する比をピークの強さ(`peak_strength`)とみなす.
7. `peak_strength` と `peak_strength_thresholds` から `quality_level` を計算する.

## `quality_level`について

`peak_strength`をs, `peak_strength_thresholds`を`[2.0, 4.0]`とした場合, 以下のように`quality_level`を算出する.

```
quality_level = 0:        s <=  2.0    Chevron Pattern確認不可
quality_level = 1:  2.0 < s <=  4.0    ぎりぎり確認可
quality_level = 2:  4.0 < s            問題なく確認可
```

## 設定ファイル

- `hampel_filter_threshold`: Hampelフィルタによるノイズ除去の閾値. この値が低いほど強くノイズが除去される.
- `hampel_filter_window_drive_frequency`: Hampelフィルタのウィンドウサイズ(入力周波数方向).
- `hampel_filter_window_time`: Hampelフィルタのウィンドウサイズ(時間方向).
- `time_spectrum_freq_band_max`: 時間方向スペクトル計算の際の, 平均を取る入力周波数方向の範囲の最大値(この値を含む).
- `peak_strength_thresholds`: quality_levelの算出に使う閾値. 詳しくは[`quality_level`について](#quality_levelについて)を参照.

## インストール

```
uv sync
```

## 実行

```
cp examples/config/config_example.json ./config.json
uv run src/main.py -c config.json -f /path/to/data.json --json
```

main.pyの出力オプション(複数可)

- `--image-dir <image_dir>`: <image_dir>に元画像・スパイクノイズ位置画像・スペクトル画像を出力する.
- `--plot`: スペクトル画像をブラウザで表示する.
- `--json`: 分析結果をjsonで出力する(詳しくは[以下](#json出力について)参照).


## json出力について

- `rabi_cycles`: Chevron Patternの時間方向の波数,
- `quality_level`: 実験結果画像の鮮明さ. 0 ~ <`peak_strength_thresholds`の要素数> の整数で, 高いほど鮮明.,
- `status`: "OK"か"ERROR". 処理中に例外が発生した場合ERROR, 単にChevron Patternが確認できなかっただけの場合はOKになる.
- `error`: 処理中に発生した例外のエラーメッセージ.

#### Chevron Patternが確認できる場合

Chevron Patternが確認できる場合, `quality_level > 0` となる.

```json
{
  "rabi_cycles": 3,
  "quality_level": 2,
  "status": "OK",
  "error": null
}
```

#### Chevron Patternが確認できない場合

Chevron Patternが確認できない場合, `quality_level = 0` となる.
ただし, 分析自体は成功しているので `status = OK` となる.
この場合でもrabi_cyclesには一応推定した値が入ることに注意.

```json
{
  "rabi_cycles": 9,
  "quality_level": 0,
  "status": "OK",
  "error": null
}
```

#### 例外発生時

コンフィグファイルの値が不正なときや, 入力ファイルにNaNが含まれている場合などに例外が発生する. `status = ERROR` となる.

```json
{
  "rabi_cycles": null,
  "quality_level": null,
  "status": "ERROR",
  "error": "zs contains NaN/Inf"
}
```

### 出力の解釈分岐

```python
if result['status'] == "ERROR":
    # 例外発生
    pass
elif result['quality_level'] == 0:
    # Chevron Patternが確認できず
    pass
else:
    # Chevron Patternが確認できた
    rabi_cycles = result['rabi_cycles']
```
