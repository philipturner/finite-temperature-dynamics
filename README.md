# Finite Temperature Dynamics

Gone are the old ways of doing quantum chemistry, where B3LYP is so computationally expensive that molecular dynamics would be unthinkable. Welcome to the era of high-quality semiempirical methods. Thanks to the same [professor](https://en.wikipedia.org/wiki/Stefan_Grimme) who invented your dispersion corrections.

Gone are the overexpensive liquid helium and cryo temperature. They're no longer needed to minimize piezo creep. We can now do mechanosynthesis at a temperature where silane doesn't freeze solid.

## Videos

[Finite temperature molecular dynamics of mechanosynthesis tooltips](https://www.youtube.com/watch?v=QEoh0x8ro_M)

[Molecular Renderer: Tutorial Video](https://www.youtube.com/watch?v=2-quQxlQWmY)

TODO: Tutorial video for reproducing the code in this repo

## Instructions

Set up Molecular Renderer:
- Install prerequisites described in [Molecular Renderer: Tutorial Video](https://www.youtube.com/watch?v=2-quQxlQWmY)
- Clone the molecular-renderer repo
- Confirm that the “Hello World” script is working
- Remove the “.git” folder from the directory
- Clone the finite-temperature-dynamics repo
- Navigate to “Sources/Workspace”
- Remove “main.swift”
- Copy over the contents of the finite-temperature-dynamics folder

Run the code:
- (macOS) Expand the Terminal window from 80 to 120 characters per line
- If monitor is 1080p instead of 4K, search for (1080, 1920) in code and replace with (540, 960)
- Run the code with the default settings
- Change background color to light blue
- Run the code again
- Change renderingOffline to true
- Change from 60 * 1 to 60 * 5 simulation frames (1 -> 5 seconds of playback)
- Run the code again
- Change tripodType and feedstockType to explore different tooltips

## Notes

Tripods shown in animation:
- HAbst
- C3Ge-H
- C3Ge-CH2
- C3Ge-SiH3
- C3Ge-CC
- C3Ge-CCBr
- N3Sn-H
- N3Sn-CH2
- N3Sn-SiH3
- N3Sn-CC
- N3Sn-CCBr

When post-processing GIF into MP4 with DaVinci Resolve, there is a major problem. The background of the video flickers, and is unbearably distracting. The problem stems from the octree quantization algorithm that compresses 2^24 possible colors to 2^8. This is the reason you see [color banding](https://en.wikipedia.org/wiki/Colour_banding) around the shiny white dot on each atom.

The default background color is grayscale, and gets binned with colors from carbon to compress the color space. This can cause flickering, as the color space changes slightly each frame. The background's assigned bin changes slightly. We solve the problem by moving to a color space region separate from any atoms in the scene. It can be remarkably close to the original color, just with a blue-enough tint to have its own spot in the color table.

In the file tree, navigate to <b>Sources</b> > <b>MolecularRenderer</b> > <b>Image</b> > <b>Render</b> > <b>RenderShader.swift</b>. Find the line of code that sets the background color to gray (0.707, 0.707, 0.707). Change this line to light blue (0.657, 0.707, 0.757). Since Molecular Renderer's `.git` folder was deleted, you won't be prompted to upstream modifications into the main repo.
