// nucleus_iod.ijm
// Nucleus detection + IOD measurement following Hardie et al. (2002)
//
// Protocol:
//   1. I_bg = 95th percentile of image
//   2. Invert (dark nuclei -> bright)
//   3. Gaussian blur sigma=4
//   4. Otsu threshold -> binary
//   5. Watershed (split touching)
//   6. Analyze Particles (size + circularity filter) - in PIXEL units
//   7. Measure IOD per nucleus on ORIGINAL (uninverted) image
//
// Run headless: imagej -b nucleus_iod.ijm "inputDir|outputFile|minSize|maxSize|minCirc|blurSigma|iBgOverride|artifactDir"
// Size args are in PIXELS (not um^2). Pipe-delimited.

// ---- PARSE ARGUMENTS ----
argString = getArgument();
if (argString == "") {
    // Fallback for Arch imagej wrapper dropping macro args in subprocess calls.
    // Python writes a temp args file before launching ImageJ.
    argsFile = "/tmp/ij_nucleus_iod_args.txt";
    if (File.exists(argsFile)) {
        argString = File.openAsString(argsFile);
        argString = replace(argString, "\n", "");
        argString = replace(argString, "\r", "");
    }
}
if (argString == "") {
    exit("Usage: imagej -b nucleus_iod.ijm 'inputDir|outputFile|minSize|maxSize|minCirc|blurSigma|iBgOverride|artifactDir'");
}

parts = split(argString, "|");
inputDir = parts[0];
outputFile = parts[1];
MIN_SIZE = parseInt(parts[2]);
MAX_SIZE = parseInt(parts[3]);
MIN_CIRC = parseFloat(parts[4]);
BLUR_SIGMA = parseFloat(parts[5]);
IBG_OVERRIDE = -1;
if (parts.length >= 7 && parts[6] != "")
    IBG_OVERRIDE = parseFloat(parts[6]);
ARTIFACT_DIR = "";
if (parts.length >= 8)
    ARTIFACT_DIR = parts[7];
MANUAL_THRESHOLD = -1;
if (parts.length >= 9) {
    if (parts[8] != "")
        MANUAL_THRESHOLD = parseInt(parts[8]);
}
PIXEL_SIZE = 0.12;
PIXEL_AREA = PIXEL_SIZE * PIXEL_SIZE;

// Ensure trailing separator
if (!endsWith(inputDir, "/")) inputDir = inputDir + "/";
if (ARTIFACT_DIR != "" && !endsWith(ARTIFACT_DIR, "/")) ARTIFACT_DIR = ARTIFACT_DIR + "/";
if (ARTIFACT_DIR != "")
    File.makeDirectory(ARTIFACT_DIR);

// ---- INIT OUTPUT ----
f = File.open(outputFile);
print(f, "filename,label,area_px,area_um2,iod,mean_od,centroid_x,centroid_y,i_bg");
File.close(f);

// ---- PROCESS ALL IMAGES ----
list = getFileList(inputDir);
imageCount = 0;
totalNuclei = 0;

setBatchMode(true);

for (idx = 0; idx < list.length; idx++) {
    fname = list[idx];
    if (!endsWith(fname, ".tiff") && !endsWith(fname, ".tif"))
        continue;

    imageCount++;
    print("Processing [" + imageCount + "] " + fname);

    open(inputDir + fname);
    originalID = getImageID();

    // NO pixel calibration - keep Analyze Particles in pixel units

    // ---- STEP 1: BACKGROUND I_bg (95th percentile) ----
    if (IBG_OVERRIDE > 0) {
        i_bg = IBG_OVERRIDE;
    } else {
        getRawStatistics(nPixels, mean, min, max, std);
        nBins = 256;
        getHistogram(values, counts, nBins);
        cumulative = 0;
        target = nPixels * 0.95;
        i_bg = max;
        for (b = 0; b < nBins; b++) {
            cumulative += counts[b];
            if (cumulative >= target) {
                i_bg = values[b];
                b = nBins;
            }
        }
    }
    print("  I_bg: " + i_bg);

    // ---- STEP 2-3: INVERT + BLUR ----
    run("Duplicate...", "title=working duplicate");
    workingID = getImageID();
    run("Invert");
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA);

    // ---- STEP 4: THRESHOLD (manual override or Otsu) ----
    if (MANUAL_THRESHOLD > 0) {
        setThreshold(MANUAL_THRESHOLD, 255);
        run("Convert to Mask");
    } else {
        setAutoThreshold("Otsu dark");
        run("Convert to Mask");
    }

    // ---- STEP 5: WATERSHED ----
    run("Watershed");

    // ---- STEP 6: ANALYZE PARTICLES (pixel units) ----
    roiManager("Reset");
    run("Analyze Particles...",
        "size=" + MIN_SIZE + "-" + MAX_SIZE +
        " circularity=" + MIN_CIRC + "-1.00" +
        " exclude add");

    nROIs = roiManager("Count");
    print("  Found " + nROIs + " nuclei");

    if (ARTIFACT_DIR != "") {
        base = replace(fname, ".tiff", "");
        base = replace(base, ".tif", "");
        selectImage(workingID);
        saveAs("Tiff", ARTIFACT_DIR + base + "__mask.tif");
        if (nROIs > 0)
            roiManager("Save", ARTIFACT_DIR + base + "__rois.zip");
    }

    // ---- STEP 7: MEASURE IOD ON ORIGINAL ----
    selectImage(originalID);

    for (r = 0; r < nROIs; r++) {
        selectImage(originalID);
        roiManager("Select", r);

        Roi.getContainedPoints(xpoints, ypoints);
        nPx = xpoints.length;
        if (nPx == 0) continue;

        // IOD = sum(log10(I_bg / I_pixel))
        iod = 0;
        sumOD = 0;
        for (p = 0; p < nPx; p++) {
            pixVal = getPixel(xpoints[p], ypoints[p]);
            if (pixVal < 1) pixVal = 1;
            od = log(i_bg / pixVal) / log(10);
            iod += od;
            sumOD += od;
        }
        meanOD = sumOD / nPx;

        // Area in um^2
        area_um2 = nPx * PIXEL_AREA;

        // Centroid from bounding box
        Roi.getBounds(rx, ry, rw, rh);
        cx = rx + rw / 2;
        cy = ry + rh / 2;

        line = fname + "," + (r + 1) + "," + nPx + "," + area_um2 + "," +
               iod + "," + meanOD + "," + cx + "," + cy + "," + i_bg;
        File.append(line, outputFile);
    }

    totalNuclei += nROIs;

    // Cleanup
    selectImage(workingID);
    close();
    selectImage(originalID);
    close();
    roiManager("Reset");
}

setBatchMode(false);
print("==== DONE ====");
print("Images: " + imageCount + ", Nuclei: " + totalNuclei);
print("Output: " + outputFile);
