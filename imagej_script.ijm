// Cleaned and documented ImageJ macro for Vimentin SunTag processing
// Purpose: Batch-process .nd2 hyperstacks for background correction,
//          generate visualizations (rainbow, gifs, tiff) and run
//          ThunderSTORM single-molecule analysis, exporting CSV results.
// Requirements:
//   - Fiji / ImageJ2 with Bio-Formats importer
//   - ThunderSTORM plugin installed
// Usage example (replace placeholders with your paths):
//   ImageJ-win64 --ij2 --console --run path/to/step1_vimentin_suntag.ijm "root_path='/path/to/images', input_folder='input', output_folder='output'"
// Notes for lay users:
//   - This script expects .nd2 files in the input folder. It will create
//     subfolders under the output folder for CSV, GIF, TIF and rainbow images.
//   - The script uses several ImageJ operations that open windows; these are
//     automatically closed when finished for each file to avoid memory leaks.

#@String root_path
#@String input_folder
#@String output_folder

// Normalize provided paths and ensure trailing separators
root_path = ensureTrailingSep(root_path);
input_folder = ensureTrailingSep(root_path + input_folder);
output_folder = ensureTrailingSep(root_path + output_folder);

// Create required output subdirectories if they do not exist
makeOutputDirs(output_folder);

// Start processing
print("Starting batch processing:\n  input_folder=" + input_folder + "\n  output_folder=" + output_folder);
process_folder(input_folder);

// Quit ImageJ when finished
run("Quit");

// ---------------------- Core processing functions ----------------------
function process_folder(folder) {
    list = getFileList(folder);
    for (i = 0; i < list.length; i++) {
        filename = list[i];
        // Skip directories and non-.nd2 files (case-insensitive)

        // Try to keep memory usage lower between files
        if (endsWith(toLowerCase(filename), ".nd2")) {
            call("java.lang.System.gc");
            print("Processing: " + folder + filename);
            process_file(folder, output_folder, filename);
        }
    }
}

function process_file(input_folder, output_folder, file) {
    // Descriptive variable names
    inputPath = input_folder + file;
    baseName = substring(file, 0, lengthOf(file) - 4); // remove .nd2

    // Paths for outputs (organized into subfolders created above)
    csvPath = output_folder + "tsoutput\\" + baseName + ".csv";
    csvsubPath = output_folder + "tsoutput\\" + baseName + "-sub.csv";
    gifPath = output_folder + "gif\\" + baseName + ".gif";
    tifPath = output_folder + "tif\\" + baseName + ".tif";
    corrGifPath = output_folder + "gif\\" + baseName + "-corr.gif";
    rainbowPath = output_folder + "rainbow\\" + baseName + ".png";

    // Open the file with Bio-Formats importer as hyperstack
    print("  Opening file: " + inputPath);
    run("Bio-Formats Importer", "open='" + inputPath + "' color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

    // Background correction: histogram matching (creates DUP_ window)
    run("Bleach Correction", "correction=[Histogram Matching]");

    // Z-projection to estimate background (creates MIN_DUP_ window)
    run("Z Project...", "projection=[Min Intensity]");

    // Subtract min-projection from corrected stack (result window created)
    imageCalculator("Subtract create stack", "DUP_" + file, "MIN_DUP_" + file);

    // The subtraction result is used for temporal color coding (rainbow)
    selectWindow("Result of DUP_" + file);
    nframes = nSlices; // number of frames in the current stack

    // Apply a temporal color code to visualize dynamics
    run("Temporal-Color Code", "lut=Spectrum start=1 end=" + nframes);

    // Save rainbow visualization and close it
    saveAs("PNG", rainbowPath);
    close();

    // Prepare corrected stack for ThunderSTORM and other outputs
    selectWindow("DUP_" + file);
    run("Duplicate...", "title=" + baseName + " duplicate range=1-" + nframes);

    // Save the corrected duplicated stack as TIFF (full precision)
    saveAs("Tiff", tifPath);

    // Enhance contrast and save corrected image as GIF preview
    run("Enhance Contrast", "saturated=0.35");
    run("8-bit");
    saveAs("Gif", corrGifPath);
    close();

    // Also save a GIF of the original (uncorrected) file for comparison
    selectWindow(file);
    run("Enhance Contrast", "saturated=0.35");
    run("8-bit");
    saveAs("Gif", gifPath);
    close();

    // Close intermediate windows used for background estimation
    selectWindow("MIN_DUP_" + file);
    close();

    // Run analysis on the subtracted image
    run("Camera setup", "offset=90.0 isemgain=false photons2adu=0.59 pixelsize=130.0");
    run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=1*std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=1.5 fitradius=3 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false renderer=[Averaged shifted histograms] magnification=1.0 colorizez=false threed=false shifts=2 repaint=5000");

    // Export ThunderSTORM results to CSV
    run("Export results", "filepath='" + csvsubPath + "' fileformat=[CSV (comma separated)] chi2=true offset=true saveprotocol=true bkgstd=true uncertainty=true intensity=true x=true sigma2=true y=true sigma1=true z=true id=true frame=true");

    // Close the average shifted histograms
    selectWindow("Averaged shifted histograms");
    close();

    // Close the subtracted image
    selectWindow("Result of DUP_" + file);
    close();

    // Run ThunderSTORM on the corrected duplicate for single-molecule localization
    // (Ensure Camera setup values match your camera; these are example values.)
    selectWindow("DUP_" + file);
    run("Camera setup", "offset=90.0 isemgain=false photons2adu=0.59 pixelsize=130.0");
    run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=1*std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=1.5 fitradius=3 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false renderer=[Averaged shifted histograms] magnification=1.0 colorizez=false threed=false shifts=2 repaint=5000");

    // Export ThunderSTORM results to CSV
    run("Export results", "filepath='" + csvPath + "' fileformat=[CSV (comma separated)] chi2=true offset=true saveprotocol=true bkgstd=true uncertainty=true intensity=true x=true sigma2=true y=true sigma1=true z=true id=true frame=true");

    // Close remaining windows for this file to free memory
    selectWindow("DUP_" + file);
    close();
    selectWindow("Averaged shifted histograms");
    close();

    print("  Finished: " + baseName);
}

// ---------------------- Small helper functions ----------------------
function ensureTrailingSep(path) {
    if (lengthOf(path) == 0)
        return path;
    if (endsWith(path, "\\") || endsWith(path, "/"))
        return path;
    return path + "\\"; // Windows-style separator (consistent with original script)
}

function makeOutputDirs(baseOutput) {
    // Create a small set of subfolders used by the script
    File.makeDirectory(baseOutput);
    File.makeDirectory(baseOutput + "tsoutput");
    File.makeDirectory(baseOutput + "gif");
    File.makeDirectory(baseOutput + "tif");
    File.makeDirectory(baseOutput + "rainbow");
}