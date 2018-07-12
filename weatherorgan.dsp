declare name "Weather Organ"; 
declare author "Mykle James Hansen"; 
declare copyright "(c) Mykle James Hansen 2018"; 
declare version "0.5"; 
declare license "This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.";

import("stdfaust.lib");

/////////////////////////
// User interface:
//
// (Note: [scale:log] sliders with min or max 0 are unuseable in some Faust UI implementations.  
// Many of the log-scale sliders below use a hack to get around this, where an integer 
// is subtracted from the slider value, to make zero values possible.
// See https://github.com/grame-cncm/faust/issues/164 for more details.
//
/////////
// Weather section:
//
// Flux (coefficient of overall flutter & drift):
flux_adj      = vslider("v:[-1]Weather/h:[0]/v:[1]/[0]Flux[unit:%]", 100, 0, 200, 1) / 100;
//
// Turbulence, in average changes per second
turbulence_adj = vslider("v:[-1]Weather/h:[0]/v:[4]/[0]Turbulence[scale:log][unit:Hz]", 0.5, 0.001, 10, 0.001);
//
// Observe flutter & drift:
flutter_dbg = vbargraph("v:[-1]Weather/h:[0]/v:[1]/[1]Flutter", -1.5, 1.5); 
drift_dbg =    vbargraph("v:[-1]Weather/h:[0]/v:[4]/[1]Drift", -1.5, 1.5); 
//
/////////
// Noise section:
//
// Noise source: White or brown noise, or gate an input signal?
noise_source_radio = hslider("v:[0]Noise/h:[0]/Source:[style:radio{'White':0;'Brown':1;'Line in':2}]", 0, 0, 2, 1);
//
// Density of noise events (inverse of "sparseness") -- from white noise to geiger counter.  Takes drift.
//density_adj 				= vslider("v:[0]Noise/h:[1]/v:[0]/[0]Density[unit:Hz][scale:log][midi:ctrl 18]", 44, 0.1, 96000, 0.1) - 0.1 : min(ma.SR); 
// NOTE: some web implementations of Faust crash when density goes below 1.0. Also, log sliders are finicky
// If compiling for the web, use this line instead:
density_adj 				= vslider("v:[0]Noise/h:[1]/v:[0]/[0]Density[unit:Hz][scale:log][midi:ctrl 18]", 44, 1, 10000, 0.1) : min(ma.SR); // web-save bounds
density_drift_adj 	= vslider("v:[0]Noise/h:[1]/v:[0]/[2]D drift[style:knob]",0,-1,1,0.001);
//
// Width: samples per noise event.  Takes flutter and drift.
width_adj 					= vslider("v:[0]Noise/h:[1]/v:[2]/[0]Width[unit:%][scale:log]", 1, 1, 101, .01) : -(1) : /(100); 
width_flutter_adj 	= vslider("v:[0]Noise/h:[1]/v:[2]/[1]W flutter[style:knob]",0,0,1,0.001);
width_drift_adj   	= vslider("v:[0]Noise/h:[1]/v:[2]/[2]W drift[style:knob]",0,-1,1,0.001);
//
// Rhythm: coefficient of periodicity for sparse_periodic_trigger.  Takes drift.
rhythm_adj 					= vslider("v:[0]Noise/h:[1]/v:[3]/[0]Rhythm[scale:exp]", 0, 0, 1, 0.01);
rhythm_drift_adj 		= vslider("v:[0]Noise/h:[1]/v:[3]/[1]R drift[style:knob]",0,-1,1,0.001);
//
// "Grit" (exponent of noise) from -1 to 1 
grit_adj 						= vslider("v:[0]Noise/h:[1]/v:[3]/[2]Grit[midi:ctrl 19]", 0, -1, 1, 0.01);
//
/////////
// Filter section
//
// Bypass filters entirely?
filter_bypass 			= checkbox("v:[1]Filter/h:[0]/Bypass");
//
// Low shelf (Beware of glitches when cutoff below 10hz)
low_shelf_adj 			= vslider("v:[1]Filter/h:[1]/[-1]Low shelf[unit:Hz][scale:log][midi:ctrl 71]", 20, 10, 22050, 1) : si.smoo;
//
// Base frequency of fundamental filter (filter 0). Takes flutter & drift.
base_center_freq 		= vslider("v:[1]Filter/h:[1]/v:[0]/[0]Freq[unit:Hz][scale:log]", 3520, 20, filter_upper_bound, 1) / 2; 	
base_flutter_adj 		= vslider("v:[1]Filter/h:[1]/v:[0]/[1]F flutter[style:knob]",0,0,1,0.001);
base_drift_adj 			= vslider("v:[1]Filter/h:[1]/v:[0]/[2]F drift[style:knob]",0,-1,1,0.001);
//
// Q of filter.  Takes flutter & drift.
Q_adj 							= vslider("v:[1]Filter/h:[1]/v:[2]/[0]Q[scale:exp][midi:ctrl 71]", 0.8, 0, 1, 0.0005);
Q_flutter_adj 			= vslider("v:[1]Filter/h:[1]/v:[2]/[1]Q flutter[style:knob]",0,0,1,0.001);
Q_drift_adj 				= vslider("v:[1]Filter/h:[1]/v:[2]/[2]Q drift[style:knob]",0,-1,1,0.001);
//
// Additional harmonic components (filters 1..N):
// Harmonic power of component (relative to base_center_freq):
filter_h_adj(N) 		= vslider("v:[1]Filter/h:[1]/h:[3]Overtones/v:[%M]/[0]O%NN ratio[style:knob]", O, 1, 15, 0.01) with {NN = N + 1; M = N + 3; O = 2*(N+1) + 1;};
// Relative amplitude of component
filter_h_level(N) 	= vslider("v:[1]Filter/h:[1]/h:[3]Overtones/v:[%M]/[1]O%NN level[style:knob]", 0, 0, 1, 0.001) with {NN = N + 1; M = N + 3;};
// Flutter and drift of component:
filter_h_flutter(N) = vslider("v:[1]Filter/h:[1]/h:[3]Overtones/v:[%M]/[2]O%NN flutter[style:knob]", 0, 0, 1, 0.001) with {NN = N + 1; M = N + 3;};
filter_h_drift_adj(N) = vslider("v:[1]Filter/h:[1]/h:[3]Overtones/v:[%M]/[3]O%NN drift[style:knob]", 0,-1,1,0.001) with {NN = N + 1; M = N + 3;};
//
/////////
// Output section:
//
// Output on/off:
outgate= checkbox("v:[99]Output/[0]Gate"); 
//
// End stage level adjustment:
outgain = vslider("v:[99]Output/[1]Gain[scale:log][midi:ctrl 7]", 1, 0.1, 10.1, 0.1) - 0.1;
//
// Optional limiter
limit_on = checkbox("v:[99]Output/[99]Limit");
//
//////////////////////////


////////////////////////
// Normalizing coefficients for noise generators:
//
// The various noise generators here and in noises.lib
// have ideal properties relative to their own theories, 
// but vary widely in absolute output level.  
normal = environment {
  no = environment { // noises.lib
    noise = 1;
    gnoise = 0.625;
    pink_noise = 12.5;
  };
  noise_pink = 0.312;
};
//
//////////////////////////


////////////////////////
// Utility functions for generation various types of sparse noise
//
//-------------`(sparse_periodic_trigger)`-----------------
// Emits +1 impulses ("trigger events") at an average frequency,
// with the distribution adjustable from purely periodic to purely random.
//
// ### Usage
// ```
// sparse_periodic_trigger(f0, periodicity, pnoise) : _
// ```
//
// Where:
//
// * `f0`: average number of triggers per second.
// * `periodicity`: coefficient of distribution noise.  0 <= periodicity <= 1.  0 = random distribution, 1 = regularly spaced pulses.
// * `pnoise`: random source of probability.  Pure white noise is good.
//----------------------------
sparse_periodic_trigger(f0, periodicity, noise) =
  (
    +(rate)               // add the rate;
    <: _, >=(1)           // if greater than 1 ...
    : _, *(1+w*noise) : - // ... subtract 1+(w*noise)
  ) ~ _
  <: _, _' : <            // emit 1 if the value decreased, 0 otherwise.
with {
  w = max(0, min(1, 1 - periodicity));
  rate = f0/ma.SR;
};
//
//
//-------------`(wide_gate)`-----------------
//-------------`(wide_hold_gate)`-----------------
// wide_gate() is a gate that opens when triggered,
// and stays open for an adjustable duration after the trigger
//
// wide_hold_gate() is similar, but instead of closing by going
// to zero, it holds the last input sample using ba.sAndH().
//
// ### Usage
//
// ```
// _ : wide_gate(trigger, width) : _
// _ : wide_hold_gate(trigger, width) : _
// ```
//
// Where:
//
// * `trigger`: The gate opens when this goes nonzero, and starts counting down when this goes to zero.
// * `width`: duration (in seconds) of the countdown before gate is closed.
//
// #### See Also
// ba.sAndH
//-------------------------------------
wide_gate(trigger, seconds, signal) = ba.if(clock > 0, signal, 0)
with {
  width = seconds * ma.SR;
  clock = max(width * (trigger != 0)) ~ (-(1) : max(0));
};
//
wide_hold_gate(trigger, seconds, signal) = ba.if(clock > 0, signal, signal:ba.sAndH(clock > 0) )
with {
  width = seconds * ma.SR;
  clock = max(width * trigger) ~ (-(1) : max(0));
};
//
//
//------------`noise_half`---------------
// Produce randomly-spaced half-impulses, at an average frequency.
// Each half-impulse is a single-sided sparse noise that doesn't return to zero;
// between events, a constant DC offset is held.
// Each event is a single change to a new amplitude.
//
// ### Usage
// ```
// noise_half(f0, pnoise, anoise) : _
// no.noise, no.gnoise(5) : noise_half(f0) : _
// no.noise <: noise_half(f0) : _
// ```
// Where:
//
// * `f0`: average number of impulses per second.
// * `pnoise`: random source of probability.  Pure white noise is good.
// * `anoise`: source of amplitude.  Pure white noise, gaussian, or any signal.
//-----------------
noise_half(f0, pnoise, anoise) = anoise : ba.sAndH(sparse_periodic_trigger(f0, 0, pnoise));
//
//
//-------------`noise_pink_sparse`-----------------
//-------------`noise_pink_sparse_n`-----------------
// Sparse "pink" noise generator.
//
// Not sure what to call this anymore, really.  The spectrum approaches pink (1/f) as
// f0 approaches the sample rate; at that point, it's using noise_half()
// to implement Trammel's "Stochastic Voss/McCartney" algorithm.
// However, at lower frequencies, the spectrum tilts toward brown (1/f^2).
//
// So it's cryptic and irregular.  As a pink noise generator it may offer some tiny
// (inadible?) amount of improved linearity over no.pink_noise, but only at great cost 
// of complexity & CPU.  It was fun and educational to write, but I have replaced it almost 
// everywhere with no.pink_noise .   However, as the perturbing force that acts on Drift, 
// it resembles (to my ear) the real behavior of wind better than anything else I've tried so far.  
// It's the right ratio of smaller movements more often and larger movements less often.  
// That's all it's used for here.
//
//
// ### Usage
//
// ```
// no.noise : noise_pink_sparse(f0) : _
// no.noise : noise_pink_sparse_n(N, f0) : _
// ```
//
// Where:
//
// * `f0`: average number of fluctuations per second.
// Set this below 20 for sparse noise (individually perceptable noise events).
// Set this to the highest value allowed by sample rate to produce general-purpose pink noise.
//
// * `N`: number of generators.  (Integer, must be known at compile time.)
// Per the Voss McCartney algorithm, full pink noise wants one generator per audible octave in the frequencies below f0.
// noise_pink_sparse defaults to 10 generators, corresponding to 10 octaves within the frequencies of 20-20000hz.
// Lower values of N may be acceptable when f0 is low, or the output is to be lowpass-filtered.
// Higher values of N, f0 and ma.SR may be appropriate for UHF listeners (bats, dolphins, audiophiles, et cetera).
//
// ### See Also:
// no.sparse_noise
// no.pink_noise_vm
//
// #### Reference:
// <http://www.ridgerat-tech.us/tech/pinkalg.htm>
//-------------------------------------
noise_pink_sparse_n(N, f0, noise) = noise, noise@N <: par(i, N,
  @(i), @(i) : noise_half(F(i))
  ) :> _
with {
  // Generator frequencies distributed exponentially between 0hz -> SR
  F(i) = ((2*f0) / (2^i)) ;     // freq decreases with i.
};
noise_pink_sparse(f0, noise) = noise_pink_sparse_n(10, f0, noise);



////////////////////////
// Sources of randomness:
//
// We use random values in various ways, and some of them sound better decorrelated.
//
// White noise sources for signal & probability should be decorrelated from each other:
s_noise = no.multinoise(2) : !, _;
p_noise = no.multinoise(2) : _, !;
//
// The drift generator and the sparse noise generator can share a random source, because they
// work at such different frequencies & do such different things ... but I still want to
// be sure they never get the same value in the same cycle.
p_drift_noise = p_noise; 	 // for triggering drift
p_sparse_noise = p_noise'; // for gating sparse noise 
//
// Gaussian noise source for sound & flutter should be decorrelated & each normalized to +-1.0:
s_gnoise = no.gnoise(5) : *(normal.no.gnoise); // sound
f_gnoise = no.gnoise(6) : *(normal.no.gnoise); // flutter values
// 
// New random values chosen at the moment of a flutter event should be decorrelated from each other:
f_gnoise_base = f_gnoise; 		// Base filter frequency
f_gnoise_base_Q = f_gnoise'; 	// Base filter Q
f_gnoise_width = f_gnoise''; 	// Noise width
f_gnoise_pan = f_gnoise''';   // Stereo position
//
// NOTE: the flutter on harmonic components are not decorrelated from f_gnoise.
// In theory I should do this, but it may change the asesthetic character of the sound,
// so it will need careful listening.
////////////////////////

////////////////////////
// Flux signals: flutter and drift
// (generated from the above noise sources.)
//
//////////
// "Flutter" is a coefficient centered around 1.0, flucutating with each noise event
// triggered by sparse_trigger.
// (It holds still between events, so the resonant filter can ring at a constant pitch.)
//
// Q_flutter: an offset centered around 0.0, +1 to -1 with a (clipped) gaussian distribution.
Q_flutter = f_gnoise_base_Q : ba.sAndH(sparse_trigger) : min(1) : max(-1) 
	: flutter_dbg //DEBUG
;
//
// width_flutter: same, but decorrelated from Q_flutter
width_flutter = f_gnoise_width : ba.sAndH(sparse_trigger) : min(1) : max(-1) ;
//
// Utility: convert flutter exponentially, so that +1 flutter raises pitch by the same interval that -1 flutter lowers it.
flutter2exponent(noise, range) = pow(2, ( noise : ba.sAndH(sparse_trigger) : *(range) ));
//
/////////
// "Drift" is low-frequency fluctuation intended to mimic the period of
// ocean waves, gusts of wind, etc.  It ranges from 0 to +1.  The rates of positive
// movement and negative movement are controlled independently, to model natural
// phenomena in which energy is delivered rapidly into a system by a particularl event,
// then leaves the system more slowly via diffusion.
// 
// The rise & decay rates of drift are a slope, not a time-bounded envelope.  Higher values of drift
// will take longer to approach and retreat.  Slope is proportionate to the overall drift rate.
// The adjusters give coefficients between 0 & 2.
// If wave goes from 0 to 1 in df seconds, the slope is (1/df) / ma.SR.  == 1 / (df * ma.SR);
// i.e. 1/samples to get there. If wave goes from 0 to 1 instantly (gets there in one sample)
// the slope is 1.0, i.e. 1/1 . The spread scross the two is 1 / (max(df * ma.SR * drift_adj, 1));
drift_slope = (2/ma.SR) * turbulence_adj;
//
// This drift signal can be thought of as a sparse brown noise wave with all negative values truncated to zero,
// passed through a lowpass-esque filter, with the rise & fall slopes determining the lowpass characteristics.
drift = flux_adj * soft_wave 
	: drift_dbg //DEBUG
with {
	// hard_wave is a sparse brown noise wave, constrained to +-1
	hard_wave = p_drift_noise : noise_pink_sparse_n(5, turbulence_adj) : *(normal.noise_pink) : min(1) : max(-1) ; 
	// soft_wave slides toward the value of hard_wave at a rate determined by drift_slope
	soft_wave = (_ <: _, slope : +) ~ _;
	slope(sig) = ba.if((sig <= hard_wave), drift_slope, -drift_slope);
};
//
// Utility: drift2exponent is like flutter2exponent, minus the sample & hold.
// coefficient from 1/2 to 2 as d_m from -1 to 1
drift2exponent(range) = pow(2, (drift * range)); 
//
/////////////////////////


/////////////////////////
// Audible noise generation:
//
// Density of noise events, with drift & flutter
// As drift varies -1/+1, vary density freq by (up to) 1/8 - 8
density = density_adj * d_drift(8 * density_drift_adj)
with {
	d_drift(range) = pow(2, (drift * range)); // coefficient from 1/2 to 2 as d_m from -1 to 1
};
//
// Adjustable sparse trigger controlled by that density.  Each trigger event sparks some sound.
sparse_trigger = sparse_periodic_trigger(density, rhythm_adj, p_sparse_noise );
//
// "Width": Adjustable gate open time (in secs) as a function of density.  Minimum 1 sample.
width = max(sw, 1/ma.SR)
with {
	// width_adj is log from 0 to 1, but scale to double-log from 0 to 1.5
	w_adj = pow(width_adj, 2) : *(1.5);
	// adjusted width: similar range.
	aw = w_adj 
				+ (width_flutter * flux_adj * width_flutter_adj / 8)
				+ (drift * flux_adj * width_drift_adj / 8)
				: max(0);
	// scaled width: divide by density
	sw = aw / density; 
};
//
// Input signal sources (choose one of 3):
//
// 1) Sparse Gaussian white noise:
gwhite_noise_source	= s_gnoise :	wide_gate(sparse_trigger, width);
//
// 2) "Sparse brown noise" (brown to pink, really; see footnote 1 in the paper for details):
brown_noise_source 	= no.pink_noise : *(normal.no.pink_noise) : wide_hold_gate(sparse_trigger, width) ;
//
// 3) Any old noise/signal from line in:
external_source 		= _ : wide_hold_gate(sparse_trigger, width);
//
// Choose one 
noise_source = gwhite_noise_source, brown_noise_source, external_source: ba.selectn(3, noise_source_radio);
//
//////////////////


/////////////////
// Timbre section -- adding character to noise events
//
///////
// "Grit": Use exponentation to push around the average levels of the random values in noise.
// a "grit" of 1.0 means noise is just noise, equally distributed & the average value would be +-.5.
// Grit is an exponent by which each sample is raised.  Numbers below 1.0 get smaller when raised to powers above 1.
// Therefore:
// exponent below 1 == compression of signal, louder output 
// exponent at 0 == max volume (all clicks are +1 or -1)
// exponent above 1 == expansion of signal, quieter output ... but this doesn't end up very useful IMO.
//
grit(noiseIn)  = abspow(noiseIn, 1 - grit_adj) 
with {
	// Samples are between [-1, 1].  We can't raise a negative number to a fractional power without imaginary numbers.
	// Instead, we raise the absolute value by the exponent, then multiply by the value's sign.
	abspow(val, power) = (abs(val) ^ power) * sign(val) ;
	sign(val) = (0 < val) - (val < 0); // -1, 0 or 1
};
//
//////////
// Fundamental (base_) Filter:
//
// The Moog VCF in vaeffects.lib creates a nice tuned resonance when the Q is set high.  
// (The three different Moog VCF implementations in vaeffects.lib  really behave 
// differently in the high notes here ... moog_vcf_2bn sounds best to me.)
resonant_filter = ve.moog_vcf_2bn;
// 
// However, moog_vcf_2bn requires some rubber bumpers:
//
// 1) There are glitch areas in this moog vcf model; some combinations of Q and freq cause obnoxious alilased noise.
// Artifacts begin to arise as F passes 22000 hz with q=0.  
// TODO: investigate: is this artifact a function of ma.SR or not?  If so, reference that.
// (Anyway, a lopass filter that starts to roll off at 22000hz is the same as a wire to human ears.)
filter_upper_bound = 22000; 
//
// 2) Negative values of Q can also produce painful dubstep, 
// and values above 1 seem to behave as compressed harmonics of the values between 0 and 1.
//
// Also, the filter bypass switch is implemented here.
filter(f, q, signal) = ba.if((filter_bypass | (f>filter_upper_bound)) , 
	signal, 
	resonant_filter(min(1, max(0, q)), f))
;
//
// Frequency of fundamental filter, with flutter & drift
base_filter_freq = base_center_freq * flutter2exponent(f_gnoise_base, flux_adj * base_flutter_adj) * drift2exponent(base_drift_adj);
//
// Q of the fundamental filter, with flutter & drift
base_filter_Q = Q_flux(exp_Q, drift, Q_drift_adj, Q_flutter, Q_flutter_adj) 
with {
	exp_Q = pow(Q_adj, 0.5); // scale is exponential, but the action is at the tip-top, so convert to double-exponential)
};
/////////////
// Overtone filters:
// 2 seems plenty, but more can easily be added here:
num_overtones = 2;
// 
// Same filter as the fundamental, with parameters computed like so:
overtone_filter(N) = filter(overtone_freq(N), overtone_Q(N)) * overtone_level(N) 
with {
	// Overtone freq = multiple of base freq + flutter & drift
	overtone_freq(N) = base_filter_freq 
		* (1+filter_h_adj(N)) 
		* overtone_f_flutter(N) 
		* overtone_drift(N); 

	// Overtone Q = same as base Q, plus flutter and drift
	overtone_Q(N) = Q_flux(base_filter_Q, drift, filter_h_drift_adj(N), Q_flutter, filter_h_flutter(N));

	// Overtone level = user control + flutter & drift
	overtone_level(N) = filter_h_level(N) 
		* overtone_drift(N)
		* overtone_f_level(N);

	// Flutter & drift for overtones:
	c_flutter(N) = Q_flutter * flux_adj * filter_h_flutter(N); // avg 0.0
	overtone_drift(N)   	= drift2exponent(flux_adj * filter_h_drift_adj(N));
	overtone_f_flutter(N) = flutter2exponent(f_gnoise, flux_adj * filter_h_flutter(N)); // TODO: DECORRELATE?
	overtone_f_level(N) 	= flutter2exponent(f_gnoise, c_flutter(N)); // TODO: DECORRELATE?
};
//
// Utility: Algorithm for applying drift and flutter to the Q of a filter, 
// while dealing with boundry conditions in a nice-sounding way.
// (TODO: use this more generally for flux on coefficients?)
//
// TODO: Clicks occur when frequency flutters too low ... since this
// filter doesn't have much effect below 40hz anyway, I'd like to cut that out.
// But rather than attenuate at the bottom, it'd be better if the flutter window
// would keep its bottom at the cutoff (40hz or whatever) and keep its set
// width for all center frequencies below (cutoff+(width/2)).  Then 
// behavior near that range would better match behavior in the rest of the
// audible range.  (I think.) (I hope.)
Q_flux(Q, drift, drift_adj, flutter_wave, flutter_adj) = Qdf : max(0) : min(1) 
with {
	drift_gap = min(Q, 1 - Q) : max(0.2); // size of the gap between the Q slider and its nearest border; minimum 20%
	// Q plus drift:
	Qd = Q + (drift * drift_adj * flux_adj * drift_gap);
	flutter_gap = min(Qd, 1 - Qd): max(0.1); // size of the flutter zone after drift is applied; minimum 10%
	// Q plus drift plus flutter:
	Qdf = Qd + (flutter_wave * flutter_adj * flux_adj * flutter_gap);
};
//
//////////////////

/////////////////
// Stereo position
//
pan_adj = vslider("v:[99]Output/[3]Pan[style:knob]", 0, -0.5, 0.5, .001);
pan_flutter_adj = vslider("v:[99]Output/[4]P flutter[style:knob]",0,0,1,0.001);
//
// pan_flutter: decorrelated from Q_flutter
pan_flutter = f_gnoise_pan : *(pan_flutter_adj) : ba.sAndH(sparse_trigger) : min(1) : max(-1) ;
pan = pan_adj + pan_flutter : min(1) : max(-1);


/////////////////
// All together now:
//
process = hgroup("[1]", 
	// Take our chosen noise source:
	noise_source 
	//
	// Delay it by 1 sample, to give the filters a head start changing parameters, to prevent a bit of crunch:
	: _'  		
	//
	// Apply grit:
	: grit 	
	//
	// Normalize levels (somewhat) between different noise sources:
	: *( (8, 1, 1) : ba.selectn(3, noise_source_radio))
	//
	// Apply low shelf :
	<: _, fi.low_shelf(-40, low_shelf_adj) : ba.if(filter_bypass) 
	//
	// Split into filter bank (fundamental plus N overtones), recombine
	<: filter(base_filter_freq, base_filter_Q), par(N, num_overtones, overtone_filter(N)) :> _
	//
	// Recombine, apply output gate
	: *(outgate) 
	//
	// Apply output gain
	: *(outgain) 
	//
	// Apply limiter
	<: _, co.limiter_1176_R4_mono : select2(limit_on) 
	// Pan stereo (normazed from +1-1 to 0-1
	: sp.panner((pan + 1)/2)

);	
