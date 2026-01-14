# Image Conversion Report

Generated: 2026-01-12 18:17:16

Source: `/run/media/jake/easystore`
Output: `data/`

## Summary Statistics

| Category | Count |
|----------|-------|
| Total output files | 292 |
| Source VSI files | 254 |
| Source BTF files | 69 |
| Valid VSI files | 252 |
| Empty VSI files (unrecoverable) | 2 |
| Valid BTF files | 63 |
| Corrupted BTF source files | 6 |
| Corrupted output files (have VSI alternative) | 2 |
| Files with disk I/O errors (recovered via dd) | 2 |
| **Unrecoverable samples** | **2** |

## Corrupted/Empty Source Files

### Empty VSI Files (0 bytes)

These files are completely empty and contain no recoverable data.

| File | Path |
|------|------|
| `ois1E0B.vsi` | `/run/media/jake/easystore` |
| `ois578B.vsi` | `/run/media/jake/easystore` |

### Corrupted BTF Files (Invalid TIFF Structure)

These files have data but corrupted TIFF headers (0 pages).

| File | Size | Path | VSI Alternative |
|------|------|------|-----------------|
| `Process_378_raw.ome.btf` | 9.7GB | `/run/media/jake/easystore/Dried_blood_nopmount/Process_378` | Yes |
| `Process_383_raw.ome.btf` | 5.9GB | `/run/media/jake/easystore/Dried_blood_nopmount/Process_383` | Yes |
| `Process_385_raw.ome.btf` | 2.6GB | `/run/media/jake/easystore/Dried_blood_nopmount/Process_385` | Yes |
| `Process_386_raw.ome.btf` | 1.9GB | `/run/media/jake/easystore/Dried_blood_nopmount/Process_386` | Yes |
| `Process_387_raw.ome.btf` | 1.2GB | `/run/media/jake/easystore/Dried_blood_nopmount/Process_387` | Yes |
| `Process_388_raw.ome.btf` | 1.8GB | `/run/media/jake/easystore/Dried_blood_nopmount/Process_388` | Yes |

### Corrupted Output Files

These output files were created but have corrupted TIFF structure (0 pages).

| Output File | Size | Source Type | Status |
|-------------|------|-------------|--------|
| `Process_342_raw_green.ome.tiff` | 2.2GB | BTF | **RECOVERED** via dd - see `recovery/` |
| `Process_344_raw_green.ome.tiff` | 2.6GB | BTF | **RECOVERED** via dd - see `recovery/` |
| `Process_360_raw_green.ome.tiff` | 622MB | BTF | Has valid VSI version (`Process_360_green.ome.tiff`) |
| `Process_361_raw_green.ome.tiff` | 491MB | BTF | Has valid VSI version (`Process_361_green.ome.tiff`) |

**Note:** The corrupted `_raw_green` files can be deleted. Process_360 and Process_361 have working VSI-converted versions.

## Unrecoverable Data Summary

The following 2 samples have no valid converted output and cannot be recovered:

| Sample | Issue | Source Files |
|--------|-------|--------------|
| `ois1E0B` | Empty VSI source (0 bytes) | `/run/media/jake/easystore/ois1E0B.vsi` |
| `ois578B` | Empty VSI source (0 bytes) | `/run/media/jake/easystore/ois578B.vsi` |

## Files with Disk I/O Errors (Unrecoverable)

The following files have valid TIFF headers but produce I/O errors when reading the image data.
Multiple conversion attempts failed at ~30% completion with `[Errno 5] Input/output error`.
This indicates bad sectors on the external drive affecting these specific files.

| Sample | Source | Source Size | Dimensions | Status |
|--------|--------|-------------|------------|--------|
| `Process_342` | `/run/media/jake/easystore/Dried_blood_nopmount/Process_342/Process_342_raw.ome.btf` | 21.8GB | 119814x59996 | **DISK ERROR** |
| `Process_344` | `/run/media/jake/easystore/Dried_blood_nopmount/Process_344/Process_344_raw.ome.btf` | 23.6GB | 94802x82417 | **DISK ERROR** |

**Recovery options:**
1. Run disk recovery tools (e.g., `ddrescue`) to attempt sector-level recovery
2. Check if backup copies exist on other media
3. Re-export from original microscope if source data is still available

## Successful Conversions

### From VSI Files

#### `/run/media/jake/easystore`

| Source File | Size | Output |
|-------------|------|--------|
| `ois1001.vsi` | 0MB | `ois1001_green.ome.tiff` |
| `ois1551.vsi` | 0MB | `ois1551_green.ome.tiff` |
| `ois19C1.vsi` | 1MB | `ois19C1_green.ome.tiff` |
| `ois1D6.vsi` | 0MB | `ois1D6_green.ome.tiff` |
| `ois1FB9.vsi` | 1MB | `ois1FB9_green.ome.tiff` |
| `ois4018.vsi` | 0MB | `ois4018_green.ome.tiff` |
| `ois4C6.vsi` | 0MB | `ois4C6_green.ome.tiff` |
| `ois5AFB.vsi` | 0MB | `ois5AFB_green.ome.tiff` |
| `ois5FFA.vsi` | 1MB | `ois5FFA_green.ome.tiff` |
| `ois601B.vsi` | 0MB | `ois601B_green.ome.tiff` |
| `ois7DEA.vsi` | 1MB | `ois7DEA_green.ome.tiff` |
| `ois98F0.vsi` | 0MB | `ois98F0_green.ome.tiff` |
| `ois9930.vsi` | 0MB | `ois9930_green.ome.tiff` |
| `oisA150.vsi` | 0MB | `oisA150_green.ome.tiff` |
| `oisA19E.vsi` | 0MB | `oisA19E_green.ome.tiff` |
| `oisA87B.vsi` | 1MB | `oisA87B_green.ome.tiff` |
| `oisB9CD.vsi` | 0MB | `oisB9CD_green.ome.tiff` |
| `oisBAB8.vsi` | 0MB | `oisBAB8_green.ome.tiff` |
| `oisBC6A.vsi` | 0MB | `oisBC6A_green.ome.tiff` |
| `oisC0AF.vsi` | 0MB | `oisC0AF_green.ome.tiff` |
| `oisD18F.vsi` | 0MB | `oisD18F_green.ome.tiff` |
| `oisD2E1.vsi` | 0MB | `oisD2E1_green.ome.tiff` |
| `oisD7F4.vsi` | 0MB | `oisD7F4_green.ome.tiff` |
| `oisD8BE.vsi` | 1MB | `oisD8BE_green.ome.tiff` |
| `oisDAC3.vsi` | 0MB | `oisDAC3_green.ome.tiff` |
| `oisEBA0.vsi` | 0MB | `oisEBA0_green.ome.tiff` |
| `oisF473.vsi` | 0MB | `oisF473_green.ome.tiff` |
| `oisF503.vsi` | 0MB | `oisF503_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1172`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_354.vsi` | 1MB | `Process_354_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1173`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_355.vsi` | 0MB | `Process_355_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1174`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_356.vsi` | 1MB | `Process_356_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1175`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_357.vsi` | 1MB | `Process_357_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1192`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_358.vsi` | 0MB | `Process_358_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1193`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_359.vsi` | 1MB | `Process_359_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1194`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_360.vsi` | 1MB | `Process_360_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1195`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_361.vsi` | 1MB | `Process_361_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1196`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_362.vsi` | 1MB | `Process_362_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1197`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_363.vsi` | 1MB | `Process_363_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1198`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_364.vsi` | 1MB | `Process_364_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1199`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_365.vsi` | 1MB | `Process_365_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1200`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_366.vsi` | 1MB | `Process_366_green.ome.tiff` |

#### `/run/media/jake/easystore/20250516/Slide glass_1203`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_367.vsi` | 1MB | `Process_367_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1207`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_368.vsi` | 1MB | `Process_368_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1208`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_369.vsi` | 0MB | `Process_369_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1209`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_370.vsi` | 1MB | `Process_370_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1210`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_371.vsi` | 0MB | `Process_371_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1211`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_372.vsi` | 0MB | `Process_372_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1212`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_373.vsi` | 1MB | `Process_373_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1213`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_374.vsi` | 0MB | `Process_374_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1214`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_375.vsi` | 0MB | `Process_375_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1215`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_376.vsi` | 1MB | `Process_376_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1216`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_377.vsi` | 1MB | `Process_377_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1217`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_378.vsi` | 0MB | `Process_378_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1218`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_379.vsi` | 0MB | `Process_379_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1219`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_380.vsi` | 1MB | `Process_380_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1220`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_381.vsi` | 0MB | `Process_381_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1221`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_382.vsi` | 1MB | `Process_382_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1222`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_383.vsi` | 0MB | `Process_383_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1224`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_384.vsi` | 1MB | `Process_384_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1228`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_385.vsi` | 1MB | `Process_385_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1229`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_386.vsi` | 1MB | `Process_386_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1230`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_387.vsi` | 1MB | `Process_387_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1249`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_388.vsi` | 0MB | `Process_388_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1250`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_389.vsi` | 1MB | `Process_389_green.ome.tiff` |

#### `/run/media/jake/easystore/20250518/Slide glass_1251`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_390.vsi` | 0MB | `Process_390_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1240`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_368.vsi` | 1MB | `Process_368_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1241`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_369.vsi` | 0MB | `Process_369_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1242`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_370.vsi` | 0MB | `Process_370_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1246`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_371.vsi` | 0MB | `Process_371_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1247`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_372.vsi` | 1MB | `Process_372_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1248`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_373.vsi` | 1MB | `Process_373_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1249`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_374.vsi` | 0MB | `Process_374_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1250`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_375.vsi` | 1MB | `Process_375_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1251`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_376.vsi` | 1MB | `Process_376_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1252`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_377.vsi` | 1MB | `Process_377_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1253`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_378.vsi` | 1MB | `Process_378_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1254`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_379.vsi` | 0MB | `Process_379_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1258`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_380.vsi` | 1MB | `Process_380_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1259`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_381.vsi` | 1MB | `Process_381_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1260`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_382.vsi` | 0MB | `Process_382_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1261`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_383.vsi` | 0MB | `Process_383_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1262`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_384.vsi` | 1MB | `Process_384_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1263`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_385.vsi` | 1MB | `Process_385_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1264`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_386.vsi` | 1MB | `Process_386_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1265`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_387.vsi` | 1MB | `Process_387_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1266`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_388.vsi` | 1MB | `Process_388_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1267`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_389.vsi` | 1MB | `Process_389_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1268`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_390.vsi` | 1MB | `Process_390_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1269`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_391.vsi` | 1MB | `Process_391_green.ome.tiff` |

#### `/run/media/jake/easystore/20250519/Slide glass_1271`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_392.vsi` | 1MB | `Process_392_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1288`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_396.vsi` | 1MB | `Process_396_green.ome.tiff` |
| `Process_397.vsi` | 0MB | `Process_397_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1289`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_393.vsi` | 0MB | `Process_393_green.ome.tiff` |
| `Process_394.vsi` | 0MB | `Process_394_green.ome.tiff` |
| `Process_395.vsi` | 0MB | `Process_395_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1290`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_398.vsi` | 1MB | `Process_398_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1291`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_399.vsi` | 0MB | `Process_399_green.ome.tiff` |
| `Process_400.vsi` | 0MB | `Process_400_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1293`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_401.vsi` | 1MB | `Process_401_green.ome.tiff` |
| `Process_402.vsi` | 1MB | `Process_402_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1294`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_403.vsi` | 1MB | `Process_403_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1295`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_404.vsi` | 1MB | `Process_404_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1296`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_405.vsi` | 0MB | `Process_405_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1297`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_406.vsi` | 0MB | `Process_406_green.ome.tiff` |
| `Process_407.vsi` | 0MB | `Process_407_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1318`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_409.vsi` | 0MB | `Process_409_green.ome.tiff` |
| `Process_410.vsi` | 0MB | `Process_410_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1319`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_411.vsi` | 1MB | `Process_411_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1320`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_412.vsi` | 0MB | `Process_412_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1321`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_413.vsi` | 1MB | `Process_413_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1322`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_414.vsi` | 0MB | `Process_414_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1323`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_415.vsi` | 1MB | `Process_415_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1324`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_416.vsi` | 1MB | `Process_416_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1325`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_417.vsi` | 1MB | `Process_417_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1326`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_418.vsi` | 1MB | `Process_418_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1327`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_419.vsi` | 1MB | `Process_419_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1329`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_420.vsi` | 0MB | `Process_420_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1330`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_421.vsi` | 1MB | `Process_421_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1331`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_422.vsi` | 1MB | `Process_422_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1338`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_423.vsi` | 1MB | `Process_423_green.ome.tiff` |
| `Process_424.vsi` | 1MB | `Process_424_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1339`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_425.vsi` | 1MB | `Process_425_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1340`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_426.vsi` | 0MB | `Process_426_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1341`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_427.vsi` | 1MB | `Process_427_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1342`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_428.vsi` | 1MB | `Process_428_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1344`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_429.vsi` | 1MB | `Process_429_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1345`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_430.vsi` | 1MB | `Process_430_green.ome.tiff` |
| `Process_431.vsi` | 1MB | `Process_431_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1347`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_432.vsi` | 1MB | `Process_432_green.ome.tiff` |
| `Process_433.vsi` | 1MB | `Process_433_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1348`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_434.vsi` | 1MB | `Process_434_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1349`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_435.vsi` | 1MB | `Process_435_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1350`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_436.vsi` | 1MB | `Process_436_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1351`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_437.vsi` | 1MB | `Process_437_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1352`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_438.vsi` | 1MB | `Process_438_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1353`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_439.vsi` | 1MB | `Process_439_green.ome.tiff` |

#### `/run/media/jake/easystore/20250521/Slide glass_1354`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_440.vsi` | 0MB | `Process_440_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1355`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_441.vsi` | 0MB | `Process_441_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1381`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_442.vsi` | 1MB | `Process_442_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1382`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_443.vsi` | 1MB | `Process_443_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1383`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_444.vsi` | 1MB | `Process_444_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1385`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_445.vsi` | 1MB | `Process_445_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1387`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_446.vsi` | 0MB | `Process_446_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1388`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_447.vsi` | 0MB | `Process_447_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1399`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_448.vsi` | 1MB | `Process_448_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1400`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_449.vsi` | 0MB | `Process_449_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1402`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_450.vsi` | 0MB | `Process_450_green.ome.tiff` |
| `Process_451.vsi` | 1MB | `Process_451_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1403`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_452.vsi` | 0MB | `Process_452_green.ome.tiff` |
| `Process_453.vsi` | 1MB | `Process_453_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1404`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_454.vsi` | 1MB | `Process_454_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1405`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_455.vsi` | 1MB | `Process_455_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1406`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_456.vsi` | 1MB | `Process_456_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1407`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_457.vsi` | 0MB | `Process_457_green.ome.tiff` |

#### `/run/media/jake/easystore/20250522/Slide glass_1408`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_458.vsi` | 1MB | `Process_458_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1409`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_459.vsi` | 1MB | `Process_459_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1410`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_460.vsi` | 0MB | `Process_460_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1411`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_461.vsi` | 0MB | `Process_461_green.ome.tiff` |
| `Process_462.vsi` | 0MB | `Process_462_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1412`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_463.vsi` | 1MB | `Process_463_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1413`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_464.vsi` | 1MB | `Process_464_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1414`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_465.vsi` | 1MB | `Process_465_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1415`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_466.vsi` | 0MB | `Process_466_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1416`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_467.vsi` | 1MB | `Process_467_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1417`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_468.vsi` | 0MB | `Process_468_green.ome.tiff` |
| `Process_469.vsi` | 0MB | `Process_469_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1418`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_470.vsi` | 1MB | `Process_470_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1419`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_471.vsi` | 1MB | `Process_471_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1420`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_472.vsi` | 0MB | `Process_472_green.ome.tiff` |
| `Process_473.vsi` | 0MB | `Process_473_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1421`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_474.vsi` | 0MB | `Process_474_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1422`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_475.vsi` | 1MB | `Process_475_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1423`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_476.vsi` | 1MB | `Process_476_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1424`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_477.vsi` | 0MB | `Process_477_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1425`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_478.vsi` | 1MB | `Process_478_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1426`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_479.vsi` | 0MB | `Process_479_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1427`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_480.vsi` | 0MB | `Process_480_green.ome.tiff` |
| `Process_481.vsi` | 1MB | `Process_481_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1428`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_482.vsi` | 1MB | `Process_482_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1429`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_483.vsi` | 0MB | `Process_483_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1430`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_484.vsi` | 0MB | `Process_484_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1431`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_485.vsi` | 0MB | `Process_485_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1435`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_486.vsi` | 1MB | `Process_486_green.ome.tiff` |
| `Process_487.vsi` | 1MB | `Process_487_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1436`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_488.vsi` | 1MB | `Process_488_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1437`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_489.vsi` | 0MB | `Process_489_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1438`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_490.vsi` | 0MB | `Process_490_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1439`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_491.vsi` | 0MB | `Process_491_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1440`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_492.vsi` | 1MB | `Process_492_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1441`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_493.vsi` | 0MB | `Process_493_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1442`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_494.vsi` | 0MB | `Process_494_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1443`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_495.vsi` | 1MB | `Process_495_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1444`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_496.vsi` | 1MB | `Process_496_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1445`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_497.vsi` | 1MB | `Process_497_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1446`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_498.vsi` | 0MB | `Process_498_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1447`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_499.vsi` | 1MB | `Process_499_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1449`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_500.vsi` | 1MB | `Process_500_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1450`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_502.vsi` | 1MB | `Process_502_green.ome.tiff` |

#### `/run/media/jake/easystore/20250523/Slide glass_1451`

| Source File | Size | Output |
|-------------|------|--------|
| `Process_501.vsi` | 1MB | `Process_501_green.ome.tiff` |

#### `/run/media/jake/easystore/_temp`

| Source File | Size | Output |
|-------------|------|--------|
| `ois12CE.vsi` | 0MB | `ois12CE_green.ome.tiff` |
| `ois1BD8.vsi` | 0MB | `ois1BD8_green.ome.tiff` |
| `ois2A40.vsi` | 0MB | `ois2A40_green.ome.tiff` |
| `ois2B27.vsi` | 0MB | `ois2B27_green.ome.tiff` |
| `ois2BA.vsi` | 0MB | `ois2BA_green.ome.tiff` |
| `ois3642.vsi` | 0MB | `ois3642_green.ome.tiff` |
| `ois3FEF.vsi` | 0MB | `ois3FEF_green.ome.tiff` |
| `ois41EF.vsi` | 0MB | `ois41EF_green.ome.tiff` |
| `ois432D.vsi` | 0MB | `ois432D_green.ome.tiff` |
| `ois526D.vsi` | 0MB | `ois526D_green.ome.tiff` |
| `ois53A5.vsi` | 1MB | `ois53A5_green.ome.tiff` |
| `ois6589.vsi` | 0MB | `ois6589_green.ome.tiff` |
| `ois6609.vsi` | 0MB | `ois6609_green.ome.tiff` |
| `ois670B.vsi` | 0MB | `ois670B_green.ome.tiff` |
| `ois6FE3.vsi` | 0MB | `ois6FE3_green.ome.tiff` |
| `ois717D.vsi` | 0MB | `ois717D_green.ome.tiff` |
| `ois751B.vsi` | 0MB | `ois751B_green.ome.tiff` |
| `ois7710.vsi` | 0MB | `ois7710_green.ome.tiff` |
| `ois79E7.vsi` | 0MB | `ois79E7_green.ome.tiff` |
| `ois8F92.vsi` | 0MB | `ois8F92_green.ome.tiff` |
| `ois955D.vsi` | 0MB | `ois955D_green.ome.tiff` |
| `ois9777.vsi` | 0MB | `ois9777_green.ome.tiff` |
| `ois98F.vsi` | 0MB | `ois98F_green.ome.tiff` |
| `ois9A5C.vsi` | 0MB | `ois9A5C_green.ome.tiff` |
| `ois9CAD.vsi` | 0MB | `ois9CAD_green.ome.tiff` |
| `ois9D5B.vsi` | 0MB | `ois9D5B_green.ome.tiff` |
| `ois9DA1.vsi` | 0MB | `ois9DA1_green.ome.tiff` |
| `ois9EE3.vsi` | 0MB | `ois9EE3_green.ome.tiff` |
| `oisA790.vsi` | 0MB | `oisA790_green.ome.tiff` |
| `oisAA88.vsi` | 0MB | `oisAA88_green.ome.tiff` |
| `oisAC8D.vsi` | 0MB | `oisAC8D_green.ome.tiff` |
| `oisAD55.vsi` | 0MB | `oisAD55_green.ome.tiff` |
| `oisAE64.vsi` | 0MB | `oisAE64_green.ome.tiff` |
| `oisAF45.vsi` | 0MB | `oisAF45_green.ome.tiff` |
| `oisB2B2.vsi` | 0MB | `oisB2B2_green.ome.tiff` |
| `oisB413.vsi` | 0MB | `oisB413_green.ome.tiff` |
| `oisB4A3.vsi` | 0MB | `oisB4A3_green.ome.tiff` |
| `oisB791.vsi` | 0MB | `oisB791_green.ome.tiff` |
| `oisBD1.vsi` | 0MB | `oisBD1_green.ome.tiff` |
| `oisC76F.vsi` | 0MB | `oisC76F_green.ome.tiff` |
| `oisCE0B.vsi` | 1MB | `oisCE0B_green.ome.tiff` |
| `oisD.vsi` | 0MB | `oisD_green.ome.tiff` |
| `oisDBB8.vsi` | 0MB | `oisDBB8_green.ome.tiff` |
| `oisDC3A.vsi` | 0MB | `oisDC3A_green.ome.tiff` |
| `oisE403.vsi` | 0MB | `oisE403_green.ome.tiff` |
| `oisE4BA.vsi` | 0MB | `oisE4BA_green.ome.tiff` |
| `oisE589.vsi` | 0MB | `oisE589_green.ome.tiff` |
| `oisE6B7.vsi` | 0MB | `oisE6B7_green.ome.tiff` |
| `oisEA29.vsi` | 0MB | `oisEA29_green.ome.tiff` |
| `oisF348.vsi` | 0MB | `oisF348_green.ome.tiff` |
| `oisF764.vsi` | 0MB | `oisF764_green.ome.tiff` |
| `oisF8B5.vsi` | 0MB | `oisF8B5_green.ome.tiff` |
| `oisFE69.vsi` | 0MB | `oisFE69_green.ome.tiff` |

### From BTF Files

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_312`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_312_raw.ome.btf` | 8.8GB | 53114x54391 | `Process_312_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_313`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_313_raw.ome.btf` | 14.5GB | 75347x63733 | `Process_313_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_314`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_314_raw.ome.btf` | 17.2GB | 78126x73075 | `Process_314_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_315`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_315_raw.ome.btf` | 10.9GB | 64230x56260 | `Process_315_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_316`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_316_raw.ome.btf` | 10.5GB | 67519x51693 | `Process_316_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_317`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_317_raw.ome.btf` | 40.2GB | 122594x108575 | `Process_317_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_318`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_318_raw.ome.btf` | 10.9GB | 78511x46040 | `Process_318_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_319`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_319_raw.ome.btf` | 13.2GB | 94902x45664 | `Process_319_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_320`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_320_raw.ome.btf` | 15.2GB | 95151x52998 | `Process_320_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_321`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_321_raw.ome.btf` | 9.3GB | 75923x40147 | `Process_321_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_322`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_322_raw.ome.btf` | 14.4GB | 102736x46392 | `Process_322_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_324`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_324_raw.ome.btf` | 17.8GB | 89243x65602 | `Process_324_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_331`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_331_raw.ome.btf` | 5.9GB | 70307x27502 | `Process_331_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_333`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_333_raw.ome.btf` | 13.8GB | 80906x56260 | `Process_333_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_334`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_334_raw.ome.btf` | 14.8GB | 80906x59996 | `Process_334_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_335`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_335_raw.ome.btf` | 15.5GB | 97581x52523 | `Process_335_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_336`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_336_raw.ome.btf` | 17.7GB | 94802x61865 | `Process_336_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_337`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_337_raw.ome.btf` | 23.6GB | 122594x63733 | `Process_337_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_338`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_338_raw.ome.btf` | 20.8GB | 114256x59996 | `Process_338_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_339`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_339_raw.ome.btf` | 10.6GB | 78126x45049 | `Process_339_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_340`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_340_raw.ome.btf` | 23.3GB | 125373x61865 | `Process_340_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_341`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_341_raw.ome.btf` | 21.8GB | 117035x61865 | `Process_341_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_342`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_342_raw.ome.btf` | 21.8GB | 119814x59996 | `Process_342_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_343`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_343_raw.ome.btf` | 16.1GB | 86464x61865 | `Process_343_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_344`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_344_raw.ome.btf` | 23.6GB | 94802x82417 | `Process_344_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_345`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_345_raw.ome.btf` | 25.2GB | 114256x73075 | `Process_345_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_346`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_346_raw.ome.btf` | 19.8GB | 80906x80549 | `Process_346_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_347`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_347_raw.ome.btf` | 7.2GB | 41997x56260 | `Process_347_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_348`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_348_raw.ome.btf` | 10.7GB | 58672x59996 | `Process_348_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_349`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_349_raw.ome.btf` | 9.2GB | 89243x33839 | `Process_349_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_350`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_350_raw.ome.btf` | 11.9GB | 61451x63733 | `Process_350_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_351`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_351_raw.ome.btf` | 20.6GB | 80906x84286 | `Process_351_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_352`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_352_raw.ome.btf` | 12.0GB | 64230x61865 | `Process_352_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_353`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_353_raw.ome.btf` | 4.7GB | 39218x39444 | `Process_353_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_354`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_354_raw.ome.btf` | 10.8GB | 61451x58128 | `Process_354_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_355`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_355_raw.ome.btf` | 12.9GB | 67010x63733 | `Process_355_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_356`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_356_raw.ome.btf` | 6.4GB | 47555x45049 | `Process_356_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_357`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_357_raw.ome.btf` | 8.4GB | 53114x52523 | `Process_357_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_358`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_358_raw.ome.btf` | 14.0GB | 67010x69338 | `Process_358_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_359`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_359_raw.ome.btf` | 10.7GB | 58672x59996 | `Process_359_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_360`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_360_raw.ome.btf` | 13.9GB | 69789x65602 | `Process_360_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_361`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_361_raw.ome.btf` | 10.8GB | 61451x58128 | `Process_361_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_362`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_362_raw.ome.btf` | 12.4GB | 64230x63733 | `Process_362_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_363`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_363_raw.ome.btf` | 17.6GB | 75347x76812 | `Process_363_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_364`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_364_raw.ome.btf` | 7.5GB | 50334x48786 | `Process_364_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_365`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_365_raw.ome.btf` | 13.6GB | 67010x67470 | `Process_365_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_366`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_366_raw.ome.btf` | 10.7GB | 58672x59996 | `Process_366_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_367`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_367_raw.ome.btf` | 9.5GB | 55893x56260 | `Process_367_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_368`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_368_raw.ome.btf` | 6.7GB | 47555x46918 | `Process_368_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_369`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_369_raw.ome.btf` | 3.3GB | 47555x22628 | `Process_369_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_370`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_370_raw.ome.btf` | 8.4GB | 53114x52523 | `Process_370_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_371`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_371_raw.ome.btf` | 3.7GB | 53114x22628 | `Process_371_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_372`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_372_raw.ome.btf` | 8.8GB | 67010x43181 | `Process_372_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_373`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_373_raw.ome.btf` | 21.0GB | 103139x67470 | `Process_373_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_374`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_374_raw.ome.btf` | 14.7GB | 69789x69338 | `Process_374_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_375`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_375_raw.ome.btf` | 14.7GB | 69789x69338 | `Process_375_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_376`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_376_raw.ome.btf` | 10.9GB | 58672x61865 | `Process_376_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_377`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_377_raw.ome.btf` | 9.5GB | 55893x56260 | `Process_377_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_379`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_379_raw.ome.btf` | 4.0GB | 36438x35707 | `Process_379_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_380`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_380_raw.ome.btf` | 5.7GB | 41997x45049 | `Process_380_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_381`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_381_raw.ome.btf` | 7.2GB | 67010x35707 | `Process_381_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_382`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_382_raw.ome.btf` | 5.5GB | 41997x43181 | `Process_382_green.ome.tiff` |

#### `/run/media/jake/easystore/Dried_blood_nopmount/Process_384`

| Source File | Size | Dimensions | Output |
|-------------|------|------------|--------|
| `Process_384_raw.ome.btf` | 6.9GB | 64230x35707 | `Process_384_green.ome.tiff` |

## BTF Files with VSI Alternatives

These BTF files were corrupted, but equivalent VSI files exist and were successfully converted.

| Process | Corrupted BTF | BTF Size | VSI Source | Output Status |
|---------|---------------|----------|------------|---------------|
| Process_378 | `Process_378/Process_378_raw.ome.btf` | 9.7GB | `Slide glass_1217/Process_378.vsi` | Converted (1.3GB) |
| Process_383 | `Process_383/Process_383_raw.ome.btf` | 5.9GB | `Slide glass_1222/Process_383.vsi` | Converted (1.1GB) |
| Process_385 | `Process_385/Process_385_raw.ome.btf` | 2.6GB | `Slide glass_1228/Process_385.vsi` | Converted (716MB) |
| Process_386 | `Process_386/Process_386_raw.ome.btf` | 1.9GB | `Slide glass_1229/Process_386.vsi` | Converted (824MB) |
| Process_387 | `Process_387/Process_387_raw.ome.btf` | 1.2GB | `Slide glass_1230/Process_387.vsi` | Converted (999MB) |
| Process_388 | `Process_388/Process_388_raw.ome.btf` | 1.8GB | `Slide glass_1249/Process_388.vsi` | Converted (1.1GB) |

