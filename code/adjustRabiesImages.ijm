function action(input, output, filename,i) {
		open(input + filename);
		run("Rotate 90 Degrees Left");
		//run("Brightness/Contrast...");
		setMinAndMax(0, 2000);
		run("Apply LUT");
		run("Size...", "width=2500 height=4399 depth=1 constrain average interpolation=Bilinear");
		outputFilename = "m1082-" + d2s(i+1,0);
		saveAs("tiff", output + filename);
		close();
}        

output = "/Volumes/GoogleDrive/My Drive/Rabies_registration/wholebrain/m1082_batch2/";
input = "/Volumes/GoogleDrive/My Drive/Confocal Images/Rabies/20200414/";

setBatchMode(true); 
list = getFileList(input);
for (i = 0; i<list.length; i++)
        action(input, output, list[i],i);
setBatchMode(false);


