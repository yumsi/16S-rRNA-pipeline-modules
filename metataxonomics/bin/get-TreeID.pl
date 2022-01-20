while (<>) {
    my $line = $_ ;
    if ($line =~ /SUCCESS\:\ (.+)/) {
	print $1 ;
    }
}
