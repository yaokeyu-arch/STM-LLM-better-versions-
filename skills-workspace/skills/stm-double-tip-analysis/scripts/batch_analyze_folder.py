#!/usr/bin/env python3
import argparse, csv, json, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_ROOT = Path(__file__).resolve().parent
ROOT = SCRIPT_ROOT.parents[2]

def get(lines, prefix):
    return next((x[len(prefix):].strip() for x in lines if x.startswith(prefix)), '')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source', type=Path)
    ap.add_argument('output', type=Path, nargs='?', help='existing/explicit run folder; default uses the saving protocol')
    ap.add_argument('--run-name', help='optional run name; otherwise use the source folder name')
    ap.add_argument('--runs-root', type=Path, default=ROOT/'runs')
    args = ap.parse_args()
    files = sorted(args.source.glob('*.sxm'))
    creator = ROOT/'skills/saving-protocol/scripts/create_run_layout.py'
    runs_root = args.output.parent if args.output is not None else args.runs_root
    requested_name = args.run_name or (args.output.name if args.output is not None else None)
    command = [sys.executable, str(creator), str(args.source), str(runs_root), '--images', *[p.name for p in files]]
    if requested_name:
        command[4:4] = ['--name', requested_name]
    created = subprocess.run(command, check=True, text=True, capture_output=True)
    run_dir = Path(created.stdout.strip())
    overall = run_dir/'overall-results'
    overall.mkdir(parents=True, exist_ok=True)
    rows=[]
    for n, source in enumerate(files, 1):
        stem=source.stem; out=run_dir/'images'/stem; out.mkdir(parents=True, exist_ok=True)
        npz=out/f'{stem}_background_subtracted.npz'
        whole=subprocess.run([sys.executable,str(SCRIPT_ROOT/'run_stm_double_tip_pipeline.py'),str(source),'-o',str(npz),'--candidates','10'],cwd=ROOT,text=True,capture_output=True)
        (out/'pipeline_stdout.log').write_text(whole.stdout)
        (out/'pipeline_stderr.log').write_text(whole.stderr)
        if whole.returncode:
            rows.append({'file':source.name,'error':whole.stderr or whole.stdout}); print(f'[{n}/{len(files)}] {source.name}: ERROR whole',flush=True); continue
        wl=whole.stdout.splitlines()
        summary_path=out/f'{stem}_sliced'/f'{stem}_sliced_summary.json'; summary=json.loads(summary_path.read_text()); scores=[float(x) for x in summary['scores']]
        top=next((x.strip() for x in wl if x.strip().startswith('(')),'')
        whole_classification=get(wl,'Classification:')
        sliced_classification=summary['classification']
        sliced_positive=sliced_classification != 'No persistent transition'
        combined_classification=('Likely double/multiple tip' if whole_classification.startswith('Likely') or sliced_positive
                                 else 'Possible/uncertain' if whole_classification.startswith('Possible')
                                 else 'No clear evidence')
        row={'file':source.name,'whole_classification':whole_classification,'whole_reason':get(wl,'Reason:'),'whole_top_candidate':top,'sliced_classification':sliced_classification,'combined_classification':combined_classification,'scan_direction':summary['scan_direction_metadata'],'acquisition_order':summary['acquisition_window_order'],'window_count':len(scores),'min_score':min(scores),'max_score':max(scores),'mean_score':sum(scores)/len(scores),'last_score':scores[-1],'change_window_index':summary['change_window_index'],'partial_correlation_window_indices':summary.get('partial_correlation_window_indices',[]),'scores':scores,'processed_npz':str(npz),'sliced_summary':str(summary_path),'window_csv':str(summary_path.with_name(f'{stem}_window_scores.csv')),'error':''}
        rows.append(row); print(f"[{n}/{len(files)}] {source.name}: {row['whole_classification']} | {row['sliced_classification']} | max={row['max_score']:.3f}",flush=True)
    (overall/'combined_analysis.json').write_text(json.dumps(rows,indent=2)+'\n')
    fields=['file','whole_classification','whole_reason','whole_top_candidate','sliced_classification','combined_classification','scan_direction','acquisition_order','window_count','min_score','max_score','mean_score','last_score','change_window_index','partial_correlation_window_indices','scores','processed_npz','sliced_summary','window_csv','error']
    with (overall/'combined_analysis.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    manifest_path = run_dir/'run_manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.update({'protocol':'saving-protocol','run_name':run_dir.name,'status':'complete','completed_utc':datetime.now(timezone.utc).isoformat(),'source_group':str(args.source.resolve()),'image_count':len(files),'successful_images':sum(not row.get('error') for row in rows),'failed_images':sum(bool(row.get('error')) for row in rows),'overall_results':['overall-results/combined_analysis.json','overall-results/combined_analysis.csv']})
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print('RUN',run_dir)
    print('DONE',overall/'combined_analysis.json',overall/'combined_analysis.csv')
if __name__=='__main__': main()
