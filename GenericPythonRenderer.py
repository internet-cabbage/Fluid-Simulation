import numpy as np
import vispy.scene
from vispy.scene import visuals, transforms
from tqdm import tqdm
from moviepy import VideoClip

filePath = "/Users/luthaisb/Code/C++/FluidDynamics/Fluid_in_box/LiquidBoxData.bin"

print("Starting to load file.")
with open(filePath, 'rb') as f:
    N, tSteps, xMax, yMax, zMax = np.fromfile(f, np.int32, count=5)
    data = np.fromfile(f, np.float32).reshape(-1, N, 3)
print("Finished loading file.")

print('Initial shape: ', data.shape)

print('N value:', N)
print('tSteps value:', tSteps)




# Make the canvas and add a simple viewer
canvas = vispy.scene.SceneCanvas(keys='interactive', show=True)
view = canvas.central_widget.add_view()

scatter = visuals.Markers(parent=view.scene)
scatter.set_data(data[0])

maxDist = np.max(np.abs(data))
#maxDist = 10
print("maxDist: ", maxDist)
view.camera = 'turntable'
view.camera.aspect = 1
view.camera.set_range(x=(-xMax, xMax), y=(-yMax, yMax), z=(-zMax, zMax))


inc = 0
stepsPerFrame = 10
def frame(event):
    global inc
    if (inc) < tSteps//stepsPerFrame:
        scatter.set_data(data[inc * stepsPerFrame])
        inc += 1
    else:
        inc = 0
        print("Looped")
print("NaN:", np.isnan(data).any(), "  Inf:", np.isinf(data).any())
timer = vispy.app.Timer(interval='auto',connect=frame,start=True)

vispy.app.run()