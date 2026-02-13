function single_flash(t; stim_start = 0.0, stim_end= 1.0, photon_flux = 1.0)
    if t >= stim_start && t <= stim_end
        return photon_flux
    else
        return 0.0
    end
end