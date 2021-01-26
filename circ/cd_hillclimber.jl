using Interpolations
using Plots
using Statistics



"""
# Secondary structure spectra standards cribbed by hand from handout
"""
wavelength = [190:1:250;]
gf_wavelength = [190, 191, 192.5, 195, 197, 200, 202, 205, 208, 210, 211,
				 214, 215, 217, 220, 222, 225, 230, 234, 238, 240, 250]	      
beta_sheet = [22400, 25300, 30000, 31900, 30000, 24300, 19300, 5700, -4700,
			  -10800, -12100, -16400, -17900, -18400, -15700, -13800, -11400,
			  -6400, 3600, -1400, 700, 0]
alpha_helix = [74800, 76900, 73300, 64300, 44300, 14300, 0, -25000., -32600,
			   -32400, -32100, -31000, -31400, -33100, -35300, -35700,
			   -32400, -21900, -11400, -4300, -3300, 0]
random_coil = [-32200.0, -34700, -37500, -41000, -41900, -36400, -25600,
			   -14500, -3400, -1400, 0, 3500, 4100, 4600, 4400, 3900, 2700,
			   800, 0, -140, -150, 0]


coil_itp = interpolate((gf_wavelength,), random_coil, Gridded(Linear()))
helix_itp = interpolate((gf_wavelength,), alpha_helix, Gridded(Linear())) 
sheet_itp = interpolate((gf_wavelength,), beta_sheet, Gridded(Linear()))

alpha_helix = helix_itp(wavelength)
beta_sheet = sheet_itp(wavelength)
random_coil = coil_itp(wavelength)

# normalize
spectral_max = maximum([abs.(alpha_helix) abs.(beta_sheet) abs.(random_coil)])
alpha_helix /= spectral_max
beta_sheet /= spectral_max
random_coil /= spectral_max

target = [14.5, 15.25, 16, 16.5, 17, 17, 17, 16.5, 16, 13, 10, 4, -2, -4, -6,
		  -8.5, -11, -11.25, -11.5, -11.625, -11.75, -11.625, -11.5, -11.25, 
		  -11, -10.75, -10.5, -10.375, -10.25, -10, -9.75, -9.5, -8.75, -8, -7,
		  -6, -5.25, -4.5, -4, -3.5, -2.75, -2, -1.75, -1.5, -1, -.5,
		  -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, -.5, 0, 0, 0, 0]

mse(aa,bb) = mean(abs.((aa-bb).^2)) 
softmax(xx) = exp.(xx) / sum(exp.(xx)) 

target = reshape(target, 61, 1)
target /= maximum(target)
components = reshape(cat(alpha_helix, beta_sheet, random_coil, dims=1), 61,3)
weights = softmax(rand(3,1))

iterations = 1000000
best_loss = Inf
losses = []
best_weights = weights
mutation_scale = 2e-1
scale_decay = 0.9995

prediction = 0

for kk in 1:iterations
	global prediction = components * weights
	global prediction /= maximum(prediction)
	loss = mse(prediction, target)
	global losses = cat(losses, loss, dims=1)

	if loss < best_loss
		global best_loss = loss
		global best_weights = weights
	end

	global weights = softmax(abs.(
			(1 - mutation_scale) * best_weights 
			+ (mutation_scale) * randn(size(weights))))
	global mutation_scale *= scale_decay

end




println("final loss after ", iterations, " steps = ", best_loss)
println("[alpha helix, beta sheet, random coil] weights: ", best_weights)

plot(prediction)
plot!(target)
savefig("./assets/prediction_target.png")
plot(losses) 
savefig("./assets/losses.png")

println("ok")
