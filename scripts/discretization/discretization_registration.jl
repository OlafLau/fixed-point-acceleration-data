# Stable integer-periodic registration; callers resample to a common grid first.
using FFTW, LinearAlgebra, Statistics

function integer_registered_metrics(a, b)
    @assert size(a)==size(b)
    aa=a.-mean(a); bb=b.-mean(b)
    fa=fft(aa); fb=fft(bb)
    peak=argmax(real.(ifft(fa.*conj.(fb))))
    shift=peak isa CartesianIndex ? Tuple(peak).-1 : (peak-1,)
    aligned=circshift(aa, .-shift)
    # Avoid cancellation in ||a||²+||b||²-2 max(correlation) near agreement.
    error=norm(aligned.-bb)/max(norm(bb),eps())
    cosine=dot(vec(abs.(fa)),vec(abs.(fb)))/max(norm(fa)*norm(fb),eps())
    return error,cosine,shift
end
