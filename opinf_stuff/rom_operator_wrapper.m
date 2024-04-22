function [Diff,Conv] = rom_operator_wrapper(options,rom_type)

options.rom.rom_type = rom_type;
options = operator_rom(options);
Diff = options.rom.Diff;
Conv = options.rom.Conv_quad;