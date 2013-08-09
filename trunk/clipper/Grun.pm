package Grun;

use Exporter;

use JSON::XS;
use File::Temp qw(tempfile);
use Carp;

our @ISA=qw(Exporter);
our @EXPORT=qw(grun grun_wait grun_kill);

my %JTMP;
sub grun_wait {
    my ($jid) = @_;
    my $ret = system("grun -q wait $jid 2>&1");
    if ($ret) {
        $ret =(($ret<<8)&255) if $ret > 255;
        $ret = 1 if !$ret;
    }

    open (my $fh, '<', $JTMP{$jid} . ".out");
    local $/=undef;
    my $out=<$fh>;
    close $fh;

    # copy pasted from below... make a function!
    $out=decode_json($out);
    if ($out->{err}) {
        # remote eval died... so we do too
        die $out->{err};
    };
    if (wantarray) {
        # return array
        return @{$out->{ret}};
    } else {
        # return single value
        return $out->{ret}->[0];
    }
}

sub grun {
    # this is only required on the execution node....so don't use it everywhere if not needed
    require B::RecDeparse;

    # at most 9 levels deep
    my $deparse=B::RecDeparse->new(level=>9);
    
    my ($op, $func, @args) = @_;
    croak("usage: grun({options}, \\\&function, \@args)") unless ref($func) eq 'CODE' && defined(wantarray);
    ($fh, $filename) = tempfile(".grun.XXXXXX", DIR=>".");
    my $code=$deparse->coderef2text($func);
    my $def=encode_json({code=>$code, args=>\@args, wantarray=>wantarray});
    print $fh $def;
    close $fh;

    my $opts;
    if ($op->{nowait}) {
        $opts = "-o $filename.out -nowait";
    }

    my $cmd = "grun $opts $^X -MGrun -e \"\\\"Grun::exec('$filename')\\\"\"";

    # get output (json string)

    my $out = `$cmd`;
    if ($op->{nowait}) {
        my ($jid) = $out =~ /job_id.*:\s*(\d+)/i;
        $JTMP{$jid}=$filename;
        return $jid;
    }

    $out=decode_json($out);
    if ($out->{err}) {
        # remote eval died... so we do too
        die $out->{err};
    };
    if (wantarray) {
        # return array
        return @{$out->{ret}};
    } else {
        # return single value
        return $out->{ret}->[0];
    }
}

sub exec {
    my ($fil) = @_;
    local $/ = undef;
    open( my $fh, '<', $fil );
    my $json = <$fh>;
    close $fh;

    my $hash=decode_json($json);
    my $sub = "sub " . $hash->{code};
    $sub = eval($sub);
    my (@ret, $ret, $err);
    eval {
        if ($hash->{wantarray}) {
            @ret=&{$sub}(@{$hash->{args}});
        } else {
            # scalar context
            $ret=&{$sub}(@{$hash->{args}});
            @ret=(($ret));
        }
    };
    my $err=$@;
    my $out=encode_json({ret=>\@ret, err=>$err});

    # return output via STDOUT
    print $out;
}
 

1;
