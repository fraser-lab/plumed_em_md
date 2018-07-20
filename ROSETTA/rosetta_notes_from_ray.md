1) Weights. There are two sampling schemes in the xml script: 1) CarteisanSampler and 2) full-atom relax in dual space (torsion and cartesian). You can think CartSampler is to sample large backbone conformational changes, and Relax is to slightly perturb the structure to drive it into an energy minimum.

My rationale behind higher density weight at Cartsampler is to make sure large backbone conformational sampling is guided by density. However, since density is not that informative around the residues we perturb, in the Relax stage (full-atom), we let Rosetta's full-atom energy function take the role to slightly perturb the loop, and reward conformations which are more protein-like.

For the comparison problem, you can simply rescore the structures with the same density weight and same scoring function.

2) If I understood correctly, you were worried that the neighboring residues might need to move to accommodate the newly sampled conformation in the loop. This part should be taken care by the energy minimizations in both CartSampler and Relax. The energy minimization is carried out for the whole protein, not just the sampled regions. Nevertheless, you can certainly play around adding movers into the protocol to see how things turn out.

idealize_jd2.linuxgccrelease -database /programs/x86_64-linux/rosetta/3.9/main/database/ -in::path ./ -in::file::s AC_SYMM.pdb -ignore_unrecognized_res -no_optH -out::path ./ -out::path::pdb ./ -chainbreaks  -overwrite

rosetta_scripts.linuxgccrelease -database /programs/x86_64-linux/rosetta/3.9/main/database/ -in::file::s AC_SYMM_ideal_edit.pdb -parser::protocol new_multi_ray.xml  -chemical::include_patches "lys_acetylated.txt" -ignore_unrecognized_res -edensity::mapreso 3.5 -default_max_cycles 200 -edensity::cryoem_scatterers -out::suffix _asymm -crystal_refine -beta
