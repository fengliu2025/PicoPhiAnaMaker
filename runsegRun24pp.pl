#!/usr/local/bin/perl
{
    $i=$ARGV[0];
    $seg=$ARGV[1];

    open(LIST, "/star/u/jjiastar/pwg/Spin/production/Run24pp_run.lst");
    @list=<LIST>;

    $run = $list[$i];
    print $run;
    chop($run);

    #system("root -b -q analy.C\\($inputlst,\\\"$rootfile\\\"\\)");
    system("root4star -b -q runPicoLambdaAnaMaker3.C\\($run,$seg\\)");
}
