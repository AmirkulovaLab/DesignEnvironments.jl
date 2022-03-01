function (pa::PressureAmplitude)(cylinders::Configuration)
    return pa(cylinders.pos)
end