import React, { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import { getAuthHeaders } from '../auth';

const STATUS_OPTIONS = [
  { label: 'All', value: '' },
  { label: 'Pending', value: 'pending' },
  { label: 'Processing', value: 'processing' },
  { label: 'Completed', value: 'completed' },
  { label: 'Failed', value: 'failed' },
];
const SORT_OPTIONS = [
  { label: 'Newest uploads', value: '-upload_time' },
  { label: 'Oldest uploads', value: 'upload_time' },
  { label: 'Latest finished', value: '-end_time' },
  { label: 'Earliest finished', value: 'end_time' },
];

function formatDate(s) {
  if (!s) return '—';
  try {
    const d = new Date(s);
    return d.toLocaleString();
  } catch {
    return s;
  }
}
function diffDuration(start, end) {
  if (!start || !end) return '—';
  const ms = new Date(end) - new Date(start);
  if (ms < 0) return '—';
  const sec = Math.floor(ms / 1000);
  const m = Math.floor(sec / 60);
  const s = sec % 60;
  return `${m}m ${s}s`;
}

function StatusTag({ status }) {
  const color =
    status === 'completed' ? '#16a34a' :
    status === 'failed' ? '#dc2626' :
    status === 'processing' ? '#2563eb' :
    '#6b7280';
  return (
    <span style={{ background: '#f3f4f6', color, border: `1px solid ${color}33`, borderRadius: 6, padding: '2px 8px', fontSize: 12 }}>
      {status}
    </span>
  );
}

function CopyId({ id }) {
  const [copied, setCopied] = useState(false);
  return (
    <button
      title="Copy Job ID"
      onClick={() => { navigator.clipboard.writeText(id); setCopied(true); setTimeout(()=>setCopied(false),1000); }}
      style={{ marginLeft: 8, fontSize: 12 }}
    >
      {copied ? '✅' : '📋'}
    </button>
  );
}

function HistoryPage() {
  async function downloadViaFetch(url, defaultName) {
    try {
      const resp = await fetch(url, { credentials: 'include', headers: getAuthHeaders() });
      if (!resp.ok) {
        const text = await resp.text();
        throw new Error(text || `HTTP ${resp.status}`);
      }
      // Try to infer filename from header; fallback to provided
      const disp = resp.headers.get('Content-Disposition') || '';
      const match = /filename=([^;]+)/i.exec(disp);
      const filename = match ? decodeURIComponent(match[1].replace(/\"/g, '')) : defaultName;
      const blob = await resp.blob();
      const link = document.createElement('a');
      link.href = URL.createObjectURL(blob);
      link.download = filename || defaultName || 'download';
      document.body.appendChild(link);
      link.click();
      setTimeout(() => {
        URL.revokeObjectURL(link.href);
        link.remove();
      }, 1000);
    } catch (e) {
      alert(`Download failed: ${e.message || e}`);
    }
  }
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);
  const [jobs, setJobs] = useState([]);
  const [total, setTotal] = useState(0);

  // Controls
  const [status, setStatus] = useState('');
  const [q, setQ] = useState('');
  const [sort, setSort] = useState('-upload_time');
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(20);

  // Drawer
  const [drawerOpen, setDrawerOpen] = useState(false);
  const [selectedJob, setSelectedJob] = useState(null);
  const [detailLoading, setDetailLoading] = useState(false);
  const [detailError, setDetailError] = useState(null);
  const [detail, setDetail] = useState(null);
  const drawerRef = useRef(null);

  function fetchList() {
    setLoading(true);
    setError(null);
    axios.get('/api/history', { params: { page, page_size: pageSize, status: status || undefined }, headers: getAuthHeaders() })
      .then(res => {
        const data = res.data || {};
        setJobs(Array.isArray(data.jobs) ? data.jobs : []);
        setTotal(Number.isFinite(data.total) ? data.total : 0);
      })
      .catch(err => {
        setError(err?.response?.data?.detail || err.message || 'Failed to load history');
      })
      .finally(() => setLoading(false));
  }

  useEffect(() => {
    fetchList();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [status, page, pageSize]);

  const filteredSorted = useMemo(() => {
    let arr = jobs;
    if (q) {
      const qq = q.toLowerCase();
      arr = arr.filter(j => (j.id?.toLowerCase().includes(qq) || j.filename?.toLowerCase().includes(qq)));
    }
    const getTime = (j, key) => (j[key] ? new Date(j[key]).getTime() : 0);
    const sorted = [...arr].sort((a,b) => {
      switch (sort) {
        case '-upload_time': return getTime(b,'upload_time') - getTime(a,'upload_time');
        case 'upload_time': return getTime(a,'upload_time') - getTime(b,'upload_time');
        case '-end_time': return getTime(b,'end_time') - getTime(a,'end_time');
        case 'end_time': return getTime(a,'end_time') - getTime(b,'end_time');
        default: return 0;
      }
    });
    return sorted;
  }, [jobs, q, sort]);

  function openDrawer(job) {
    setSelectedJob(job);
    setDrawerOpen(true);
    setDetail(null);
    setDetailError(null);
    setDetailLoading(true);
    axios.get(`/api/history/${job.id}`, { headers: getAuthHeaders() })
      .then(async res => {
        const d = res.data;
        // If has umap plots, try to build thumbnails by reading the UMAP ZIP
        const hasUmap = (d.artifacts || []).some(a => a.a_type === 'umap_plot');
        if (hasUmap) {
          try {
            const JSZip = (await import('jszip')).default;
            const zipResp = await axios.get(`/api/download/${job.id}/umap`, { responseType: 'arraybuffer', headers: getAuthHeaders() });
            const zip = await JSZip.loadAsync(zipResp.data);
            const pngFiles = Object.keys(zip.files).filter(name => name.endsWith('.png'));
            const thumbs = [];
            for (const name of pngFiles.slice(0, 12)) {
              const file = zip.files[name];
              const blob = await file.async('blob');
              const url = URL.createObjectURL(blob);
              thumbs.push({ name, url, blob });
            }
            d.__thumbs = thumbs;
          } catch (e) {
            // ignore thumbnail errors, keep other details
            d.__thumbs = [];
          }
        }
        setDetail(d);
      })
      .catch(err => setDetailError(err?.response?.data?.detail || err.message || 'Failed to load job detail'))
      .finally(() => setDetailLoading(false));
  }

  useEffect(() => {
    function onKey(e) {
      if (e.key === 'Escape') setDrawerOpen(false);
    }
    if (drawerOpen) {
      document.addEventListener('keydown', onKey);
      // focus trap
      setTimeout(() => drawerRef.current?.focus(), 0);
    } else {
      document.removeEventListener('keydown', onKey);
    }
    return () => document.removeEventListener('keydown', onKey);
  }, [drawerOpen]);

  const pageCount = Math.max(1, Math.ceil(total / pageSize));

  return (
    <div style={{ padding: '16px' }}>
      <div className="card">
        {/* Controls */}
        <div style={{ display: 'flex', gap: 12, alignItems: 'center', marginBottom: 12, flexWrap: 'wrap' }}>
          <label> Status:
            <select value={status} onChange={e => { setPage(1); setStatus(e.target.value); }} style={ctl}>
              {STATUS_OPTIONS.map(o => <option key={o.value} value={o.value}>{o.label}</option>)}
            </select>
          </label>
          <label> Search:
            <input value={q} onChange={e=>setQ(e.target.value)} placeholder="Job ID or filename" style={ctl} />
          </label>
          <label> Sort:
            <select value={sort} onChange={e=>setSort(e.target.value)} style={ctl}>
              {SORT_OPTIONS.map(o => <option key={o.value} value={o.value}>{o.label}</option>)}
            </select>
          </label>
          <label> Page size:
            <select value={pageSize} onChange={e=>{ setPageSize(Number(e.target.value)); setPage(1); }} style={ctl}>
              {[10,20,50].map(n=> <option key={n} value={n}>{n}</option>)}
            </select>
          </label>
          <div style={{ marginLeft: 'auto' }}>
            <button onClick={fetchList}>↻ Refresh</button>
          </div>
        </div>

        {/* Error banner */}
        {error && (
          <div style={{ background:'#fee2e2', color:'#b91c1c', padding:8, border:'1px solid #fecaca', borderRadius:6, marginBottom:12 }}>
            <strong>Error:</strong> {String(error)} <button onClick={fetchList} style={{ marginLeft: 8 }}>Retry</button>
          </div>
        )}

        {/* Table */}
        {loading ? (
          <div className="skeleton-table">Loading…</div>
        ) : filteredSorted.length === 0 ? (
          <div style={{ opacity: 0.8 }}>No results. Try changing filters or search.</div>
        ) : (
          <div style={{ overflowX: 'auto' }}>
            <table style={{ width: '100%', borderCollapse: 'collapse', background: 'white' }}>
              <thead>
                <tr>
                  <th style={thStyle}>Job ID</th>
                  <th style={thStyle}>Dataset</th>
                  <th style={thStyle}>Status</th>
                  <th style={thStyle}>Upload</th>
                  <th style={thStyle}>End</th>
                  <th style={thStyle}>Duration</th>
                  <th style={thStyle}>Actions</th>
                </tr>
              </thead>
              <tbody>
                {filteredSorted.map(job => {
                  const jobId = job.id;
                  const csvUrl = `/api/download/${jobId}/file`;
                  const umapUrl = `/api/download/${jobId}/umap`;
                  const datasetName = job.dataset?.name || job.filename || '—';
                  return (
                    <tr key={jobId} onClick={() => openDrawer(job)} style={{ cursor: 'pointer' }}>
                      <td style={tdStyle}><code>{jobId}</code><CopyId id={jobId} /></td>
                      <td style={tdStyle}>{datasetName}</td>
                      <td style={tdStyle}><StatusTag status={job.status} /></td>
                      <td style={tdStyle}>{formatDate(job.upload_time)}</td>
                      <td style={tdStyle}>{formatDate(job.end_time)}</td>
                      <td style={tdStyle}>{diffDuration(job.start_time, job.end_time)}</td>
                      <td style={tdStyle} onClick={(e)=>e.stopPropagation()}>
                        {job.csv_available ? (
                          <button
                            onClick={(e)=>{ e.stopPropagation(); downloadViaFetch(csvUrl, `${jobId}.csv`); }}
                            style={{ marginRight: 8 }}
                          >
                            CSV
                          </button>
                        ) : (
                          <button disabled style={{ marginRight: 8 }}>CSV</button>
                        )}
                        {job.umap_available ? (
                          <button
                            onClick={(e)=>{ e.stopPropagation(); downloadViaFetch(umapUrl, `umap_results_${jobId}.zip`); }}
                          >
                            UMAP ZIP
                          </button>
                        ) : (
                          <button disabled>UMAP ZIP</button>
                        )}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        )}

        {/* Pagination */}
        <div style={{ display:'flex', alignItems:'center', gap:12, marginTop:12 }}>
          <button disabled={page<=1} onClick={()=>setPage(p=>Math.max(1,p-1))}>Prev</button>
          <span>Page {page} / {pageCount}</span>
          <button disabled={page>=pageCount} onClick={()=>setPage(p=>Math.min(pageCount,p+1))}>Next</button>
          <span style={{ marginLeft: 'auto' }}>Total: {total}</span>
        </div>
      </div>

      {/* Drawer */}
      {drawerOpen && (
        <div className="drawer-backdrop" onClick={()=>setDrawerOpen(false)}>
          <div
            className="drawer"
            role="dialog"
            aria-modal="true"
            tabIndex={-1}
            ref={drawerRef}
            onClick={(e)=>e.stopPropagation()}
          >
            <div style={{ display:'flex', justifyContent:'space-between', alignItems:'center' }}>
              <h3 style={{ margin:0 }}>Job {selectedJob?.id}</h3>
              <button onClick={()=>setDrawerOpen(false)} aria-label="Close">✖</button>
            </div>
            {detailLoading ? (
              <div className="skeleton">Loading detail…</div>
            ) : detailError ? (
              <div style={{ background:'#fee2e2', color:'#b91c1c', padding:8, border:'1px solid #fecaca', borderRadius:6 }}>
                {String(detailError)}
              </div>
            ) : detail ? (
              <div style={{ marginTop: 8 }}>
                <div style={{ marginBottom: 8 }}>
                  <div><strong>Status:</strong> <StatusTag status={detail.status} /></div>
                  <div><strong>Upload:</strong> {formatDate(detail.upload_time)}</div>
                  <div><strong>Start:</strong> {formatDate(detail.start_time)}</div>
                  <div><strong>End:</strong> {formatDate(detail.end_time)}</div>
                  <div><strong>Duration:</strong> {diffDuration(detail.start_time, detail.end_time)}</div>
                  {detail.dataset?.name && <div><strong>Dataset:</strong> {detail.dataset.name}</div>}
                </div>
                <div style={{ marginBottom: 8 }}>
                  <button
                    onClick={()=>downloadViaFetch(`/api/download/${detail.id}/file`, `${detail.id}.csv`)}
                    style={{ marginRight: 8 }}
                  >
                    Download CSV
                  </button>
                  <button
                    onClick={()=>downloadViaFetch(`/api/download/${detail.id}/umap`, `umap_results_${detail.id}.zip`)}
                  >
                    Download UMAP ZIP
                  </button>
                </div>
                {detail.__thumbs && detail.__thumbs.length > 0 && (
                  <div>
                    <h4>UMAP Plots</h4>
                    <div style={{ display:'grid', gridTemplateColumns:'repeat(auto-fill, minmax(120px, 1fr))', gap:8 }}>
                      {detail.__thumbs.map((t,i)=> (
                        <a key={i} href={t.url} target="_blank" rel="noreferrer">
                          <img src={t.url} alt={t.name} style={{ width:'100%', height:100, objectFit:'cover', border:'1px solid #eee', borderRadius:6 }} loading="lazy" />
                        </a>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            ) : null}
          </div>
        </div>
      )}
    </div>
  );
}

const ctl = { marginLeft: 6 };
const thStyle = { textAlign: 'left', borderBottom: '1px solid #ddd', padding: '8px' };
const tdStyle = { borderBottom: '1px solid #f0f0f0', padding: '8px', verticalAlign: 'top' };

export default HistoryPage;


