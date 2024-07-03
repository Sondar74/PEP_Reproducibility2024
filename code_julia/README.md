# Benjamin-Bona-Mahony (BBM) Equation

The code in this directory is used to create the data for the BBM equation.
First, start Julia in this directory. On Linux and macOS, you can do this
by running

```bash
julia
```

in the terminal. Then, you can execute the following code in the Julia REPL:

```julia
julia> include("code_bbm.jl")

julia> main()
```

This will print some output to the screen.




1.  Launch a Julia prompt and paste/execute the contents of `code_bbm.jl`.
This will generate arrays containing the data used for the plot.

Then execute:

    bbm_plot(false)

to generate the plot of the Explicit-Midpoint method, and execute

    bbm_plot(true)

to generate the plot of the Pseudo-Energy-Preserving method "PEP(4,2,5)".

2. Execute

    generate_arrays(false)

    generate_arrays(true)

to save the arrays containing the data used for the plots of the RK(2,2) and PEP(4,2,5) methods respectively.

3. Open `BBM_plots.ipynb` and execute all cells of the notebook.

This code has been tested using Julia version 1.6.0.
