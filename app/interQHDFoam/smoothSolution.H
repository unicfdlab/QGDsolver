if (args.found("smoothAlpha"))
{
    const label nSmoothIters = args.getOrDefault<label>("nSmoothIters", 1);
    const scalar smoothCoeff = args.getOrDefault<scalar>("smoothCoeff", 0.75);
    
    for (label nSmoothIter = 0; nSmoothIter < nSmoothIters; nSmoothIter++)
    {
        Info << "Smooth iteration: " << nSmoothIter << endl;
        fvc::smooth(alpha1, smoothCoeff);
    }
}