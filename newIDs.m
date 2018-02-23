function IDs = newIDs()
% 0 is currently prepended just to conform to the iteration standard of the other implementation
IDs = cell(2^9 * 2^3,1);
%disp(length(IDs))
idx = 1;
for Q = 1:2
	for AB = 1:2
		for c = 0:2
			for gr = 0:2
				for gi = 0:1
					% cross u real terms
					for Omr = 0:1
						% imaginary theta terms
						for Gami = 0:1
							% real K term
							for Kr = 0:1
								% cross u imag terms
								for Omi = 0:1
									% real theta terms
									for Gamr = 0:1
										% imag K term
										for Ki = 0:1
											%%% exceptions to toss out
											if (AB == 2) && ((Omi == 1) || (Gamr == 1))
												% no change for real parts of Gamma and imaginary of u for simplification
												continue
											elseif (gi == 1) && (Omi == 1)
												% no gs in imaginary part of simplified u
												continue
											elseif (AB == 2) && (gi == 1) && (Omi == 0)
												% no gs ever in imaginary part of u under option b
												continue
											elseif (AB==2) && (gr == 1) && (Omr == 1)
												% Omega_r is 0 when gr = 1 under option 2
												continue
											end
											IDs{idx} = sprintf('0%d%d%d%d%d%d%d%d%d%d%d',...
																Q, AB, c, gr, gi, Omr, Gami, Kr, Omi, Gamr, Ki );
											idx = idx + 1;
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end
end
IDs = IDs(1:idx-1);
